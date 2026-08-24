using System.Threading;

namespace SkylineCadenza.Core.Scheduling;

/// <summary>
/// Greedy load-balanced MTM target selection. Direct port of
/// <c>schedule_targets</c> from
/// <c>targeted-modeling/gpf_coverage.ipynb</c> (cell <c>cc94d293</c>).
/// </summary>
/// <remarks>
/// Two passes:
/// <list type="number">
/// <item>
/// <b>Cover</b>: process protein groups in order of ascending candidate
/// count; for each, attempt the best peptide, drop it if it doesn't fit,
/// try the next, etc., until one succeeds or the group is exhausted.
/// Smallest groups are tried first because they have no fallbacks.
/// </item>
/// <item>
/// <b>Load up</b>: round-robin over already-covered groups, attempting to
/// add the next-best peptide. Loop until no group adds anything.
/// </item>
/// </list>
///
/// A peptide fits if it can either join an existing slot (m/z range stays
/// within window cap, RT co-elution intersection stays non-empty, top-4
/// fragments don't clash, and any RT extension keeps every newly-covered
/// bin within budget) or open a new slot inside budget. The intersection-
/// based co-elution requirement enforces "all multiplexed precursors are
/// genuinely sampled at the same instant" rather than just having
/// touching bounding boxes.
/// </remarks>
public static class Scheduler
{
    public static ScheduleResult Run(
        IReadOnlyList<Candidate> candidates,
        SchedulingParameters parameters,
        CancellationToken cancellationToken = default)
    {
        if (candidates.Count == 0)
        {
            return new ScheduleResult
            {
                ScheduledIndices = Array.Empty<int>(),
                ScheduledSlotIds = Array.Empty<int>(),
                Slots = Array.Empty<Slot>(),
                RtGrid = Array.Empty<double>(),
                SlotCountCurve = Array.Empty<int>(),
                ProteinGroupsCovered = 0,
            };
        }

        double padMin = parameters.FiringPadMin;
        // Grid origin and bin alignment are derived from the EXPERIMENTAL
        // candidates only, so the experimental tier's per-bin budget checks -
        // and therefore the whole experimental schedule - are bit-identical
        // whether or not a predicted (backfill) tier is present. rtHi is then
        // extended over ALL candidates so predicted peptides at later RT still
        // land in a valid bin; predicted peptides earlier than the first
        // experimental one clamp into bin 0 (RtToBinRange floors at 0), which is
        // acceptable because backfill rarely precedes the measured gradient. A
        // pool with no experimental candidates falls back to the full range.
        double rtLo = double.PositiveInfinity;
        double rtHi = double.NegativeInfinity;
        bool anyExperimental = false;
        foreach (var c in candidates)
        {
            if (c.Source != CandidateSource.Experimental) continue;
            anyExperimental = true;
            if (c.RtStart < rtLo) rtLo = c.RtStart;
            if (c.RtStop > rtHi) rtHi = c.RtStop;
        }
        if (!anyExperimental)
        {
            foreach (var c in candidates)
            {
                if (c.RtStart < rtLo) rtLo = c.RtStart;
                if (c.RtStop > rtHi) rtHi = c.RtStop;
            }
        }
        // Extend the upper bound over predicted candidates without moving the
        // experimental-derived origin.
        foreach (var c in candidates)
            if (c.RtStop > rtHi) rtHi = c.RtStop;
        rtLo = rtLo - padMin - 0.1;
        rtHi = rtHi + padMin + 0.1;
        int nBins = (int)Math.Ceiling((rtHi - rtLo) / parameters.RtBinMin);
        var slotsPerBin = new int[nBins];

        (int a, int b) RtToBinRange(double rtStart, double rtStop)
        {
            int a = (int)((rtStart - rtLo) / parameters.RtBinMin);
            if (a < 0) a = 0;
            int b = (int)((rtStop - rtLo) / parameters.RtBinMin) + 1;
            if (b > nBins) b = nBins;
            return (a, b);
        }

        int RangeMax(int[] arr, int a, int b)
        {
            int max = int.MinValue;
            for (int i = a; i < b; i++)
                if (arr[i] > max) max = arr[i];
            return max;
        }

        void RangeAddOne(int[] arr, int a, int b)
        {
            for (int i = a; i < b; i++) arr[i]++;
        }

        var slots = new List<Slot>();
        var slotByMzBin = new Dictionary<int, List<int>>();

        void RegisterSlotMzBins(int sid, double mzMin, double mzMax)
        {
            int lo = (int)Math.Floor(mzMin);
            int hi = (int)Math.Floor(mzMax);
            for (int b = lo; b <= hi; b++)
            {
                if (!slotByMzBin.TryGetValue(b, out var list))
                {
                    list = new List<int>();
                    slotByMzBin[b] = list;
                }
                if (!list.Contains(sid)) list.Add(sid);
            }
        }

        // Slot-aware packing (WindowPacking.Disjoint). Two guards, with
        // different lifetimes:
        //
        //   1. Growth bounds from CLASH BARRIERS. Sourced from the candidate
        //      POOL, not from the slots built so far - a precursor co-isolates
        //      whenever it falls inside a firing band, whether or not it ever
        //      wins a slot. Because the pool is static, the bound a slot gets
        //      is independent of the order proteins happen to be processed in.
        //   2. NEIGHBOUR BANDS. Slots appear as scheduling proceeds, so these
        //      are probed at expansion time rather than cached.
        bool disjointPacking = parameters.Packing == WindowPacking.Disjoint
                               && parameters.Mode == AcquisitionMode.Mtm;

        // Clash-barrier pool for TightenBounds. Populated per-tier inside the
        // scheduling loop below: the experimental tier sees only experimental
        // barriers (so its slot growth bounds don't depend on the predicted
        // pool), the predicted tier sees all of them. Empty until then.
        var candByMzBin = new Dictionary<int, List<int>>();

        double prmHalf = parameters.PrmIsolationWidthTh / 2.0;
        // Furthest a barrier can sit and still be swallowed: the widest the
        // member centers may span, plus the half-quad the band reaches past
        // the edge member. Anything beyond this is unreachable, so the scan
        // stays bounded.
        double barrierReach =
            Math.Max(0.0, parameters.IsolationWindowTh - parameters.PrmIsolationWidthTh) + prmHalf;

        // Tighten a slot's growth bounds with the barriers one member
        // introduces: the nearest co-eluting precursor on each side whose
        // fragments clash with it. The emitted band reaches prmHalf past the
        // edge member, so the bound sits half a quad width inside the barrier.
        // Bounds only ever tighten, so re-running this as members are added is
        // monotone and order-independent.
        void TightenBounds(Slot slot, int memberIdx)
        {
            if (!disjointPacking) return;
            var m = candidates[memberIdx];
            int lo = (int)Math.Floor(m.PrecursorMz - barrierReach);
            int hi = (int)Math.Floor(m.PrecursorMz + barrierReach);
            for (int b = lo; b <= hi; b++)
            {
                if (!candByMzBin.TryGetValue(b, out var idxs)) continue;
                foreach (int xi in idxs)
                {
                    if (xi == memberIdx || slot.MemberIndices.Contains(xi)) continue;
                    var x = candidates[xi];
                    if (Math.Abs(x.PrecursorMz - m.PrecursorMz) > barrierReach) continue;
                    // Only a precursor that co-elutes with this member can be
                    // co-isolated with it.
                    if (!(x.RtStart < m.RtStop && m.RtStart < x.RtStop)) continue;
                    if (!FragmentClash.AnyWithin(x.Top4Fragments, m.Top4Fragments, parameters.FragmentTolDa))
                        continue;

                    if (x.PrecursorMz > slot.MzMax)
                        slot.GrowHi = Math.Min(slot.GrowHi, x.PrecursorMz - prmHalf);
                    else if (x.PrecursorMz < slot.MzMin)
                        slot.GrowLo = Math.Max(slot.GrowLo, x.PrecursorMz + prmHalf);
                    // A barrier already inside [MzMin, MzMax] cannot be excluded
                    // by any band around these members; nothing to bound.
                }
            }
        }

        // True when a band firing over [rtStart, rtStop] would overlap the band
        // of a DIFFERENT, concurrently-active slot. Pass selfSid = -1 when
        // testing a slot that does not exist yet.
        bool OverlapsOtherSlotBand(int selfSid, double bandLo, double bandHi, double rtStart, double rtStop)
        {
            if (!disjointPacking) return false;
            int lo = (int)Math.Floor(bandLo) - 1;
            int hi = (int)Math.Floor(bandHi) + 1;
            for (int b = lo; b <= hi; b++)
            {
                if (!slotByMzBin.TryGetValue(b, out var sids)) continue;
                foreach (int sid in sids)
                {
                    if (sid == selfSid) continue;
                    var s = slots[sid];
                    double sLo = s.MzMin - prmHalf, sHi = s.MzMax + prmHalf;
                    if (bandLo < sHi && sLo < bandHi && rtStart < s.RtStop && s.RtStart < rtStop)
                        return true;
                }
            }
            return false;
        }

        // Read-only feasibility check: true iff `candIndex` could be merged
        // into some existing MTM slot right now without opening a new one.
        // Mirrors the merge-half of TrySchedule and MUST be kept in sync
        // with it. Used by the MaximizeProteins cover pass to prefer
        // peptides that don't consume a slot.
        bool CanJoinExistingSlot(int candIndex)
        {
            if (parameters.Mode == AcquisitionMode.Prm) return false;
            var prec = candidates[candIndex];
            var frags = prec.Top4Fragments;
            if (frags.Length == 0) return false;

            double mz = prec.PrecursorMz;
            double rtStart = prec.RtStart;
            double rtStop = prec.RtStop;
            double rtStartPad = rtStart - padMin;
            double rtStopPad = rtStop + padMin;
            double windowTh = parameters.IsolationWindowTh;
            // Slot-edge rule: see TrySchedule for the derivation. Every
            // member's PrmIsolationWidthTh quadrupole window must fit
            // inside the slot, so center spans are capped at
            // (windowTh - PrmIsolationWidthTh).
            double centerSpanBudget = Math.Max(0.0, windowTh - parameters.PrmIsolationWidthTh);

            int mzBinLo = (int)Math.Floor(mz - windowTh);
            int mzBinHi = (int)Math.Floor(mz + windowTh);
            var seen = new HashSet<int>();
            for (int bin = mzBinLo; bin <= mzBinHi; bin++)
            {
                if (!slotByMzBin.TryGetValue(bin, out var sids)) continue;
                foreach (int sid in sids)
                {
                    if (!seen.Add(sid)) continue;
                    var slot = slots[sid];

                    double newMzMin = Math.Min(slot.MzMin, mz);
                    double newMzMax = Math.Max(slot.MzMax, mz);
                    if (newMzMax - newMzMin > centerSpanBudget) continue;
                    // Growth bounds: never expand across a clashing precursor.
                    if (newMzMin < slot.GrowLo || newMzMax > slot.GrowHi) continue;

                    double newCoStart = Math.Max(slot.CoStart, rtStart);
                    double newCoStop = Math.Min(slot.CoStop, rtStop);
                    if (newCoStart >= newCoStop) continue;

                    if (FragmentClash.AnyWithin(frags, slot.Fragments, parameters.FragmentTolDa))
                        continue;

                    if (parameters.ChargeHandling == ChargeHandling.SameChargePerSlot
                        && slot.MemberIndices.Count > 0
                        && candidates[slot.MemberIndices[0]].PrecursorCharge != prec.PrecursorCharge)
                        continue;

                    double newRtStartPad = Math.Min(slot.RtStart, rtStartPad);
                    double newRtStopPad = Math.Max(slot.RtStop, rtStopPad);
                    var (oldA, oldB) = RtToBinRange(slot.RtStart, slot.RtStop);
                    var (extA, extB) = RtToBinRange(newRtStartPad, newRtStopPad);
                    int extLoFrom = extA, extLoTo = Math.Min(oldA, extB);
                    int extHiFrom = Math.Max(oldB, extA), extHiTo = extB;
                    if (extLoFrom < extLoTo && RangeMax(slotsPerBin, extLoFrom, extLoTo) + 1 > parameters.CycleBudget)
                        continue;
                    if (extHiFrom < extHiTo && RangeMax(slotsPerBin, extHiFrom, extHiTo) + 1 > parameters.CycleBudget)
                        continue;
                    if (OverlapsOtherSlotBand(sid, newMzMin - prmHalf, newMzMax + prmHalf,
                                              newRtStartPad, newRtStopPad))
                        continue;

                    return true;
                }
            }
            return false;
        }

        // Try to schedule a single precursor. Returns the slot id it was
        // assigned to, or -1 if it didn't fit.
        int TrySchedule(int candIndex)
        {
            var prec = candidates[candIndex];
            double mz = prec.PrecursorMz;
            double rtStart = prec.RtStart;
            double rtStop = prec.RtStop;
            double rtStartPad = rtStart - padMin;
            double rtStopPad = rtStop + padMin;
            var frags = prec.Top4Fragments;
            double windowTh = parameters.Mode == AcquisitionMode.Prm ? 0.0 : parameters.IsolationWindowTh;
            // Slot-edge rule: every member's quadrupole window
            // (PrmIsolationWidthTh wide, centered on the precursor m/z)
            // must fit fully inside the slot's nominal isolation window.
            // Equivalently, every member's center must sit at least
            // PrmIsolationWidthTh / 2 from each edge, so the max span of
            // precursor centers in one slot is (windowTh -
            // PrmIsolationWidthTh).
            double centerSpanBudget = Math.Max(0.0, windowTh - parameters.PrmIsolationWidthTh);

            if (parameters.Mode == AcquisitionMode.Mtm && frags.Length > 0)
            {
                int mzBinLo = (int)Math.Floor(mz - windowTh);
                int mzBinHi = (int)Math.Floor(mz + windowTh);
                var seen = new HashSet<int>();
                for (int bin = mzBinLo; bin <= mzBinHi; bin++)
                {
                    if (!slotByMzBin.TryGetValue(bin, out var sids)) continue;
                    foreach (int sid in sids)
                    {
                        if (!seen.Add(sid)) continue;
                        var slot = slots[sid];

                        double newMzMin = Math.Min(slot.MzMin, mz);
                        double newMzMax = Math.Max(slot.MzMax, mz);
                        if (newMzMax - newMzMin > centerSpanBudget) continue;
                        // Growth bounds: never expand across a clashing precursor.
                        if (newMzMin < slot.GrowLo || newMzMax > slot.GrowHi) continue;

                        // Strict co-elution: intersection of unpadded peak
                        // ranges must remain non-empty.
                        double newCoStart = Math.Max(slot.CoStart, rtStart);
                        double newCoStop = Math.Min(slot.CoStop, rtStop);
                        if (newCoStart >= newCoStop) continue;

                        if (FragmentClash.AnyWithin(frags, slot.Fragments, parameters.FragmentTolDa))
                            continue;

                        // Same-charge constraint: skip slots whose existing
                        // members carry a different precursor charge. With
                        // SameChargePerSlot every member shares the first
                        // member's charge, so one comparison is enough.
                        if (parameters.ChargeHandling == ChargeHandling.SameChargePerSlot
                            && slot.MemberIndices.Count > 0
                            && candidates[slot.MemberIndices[0]].PrecursorCharge != prec.PrecursorCharge)
                            continue;

                        // Extending the padded firing window. The set of bins
                        // already incremented for this slot is the half-open
                        // range [oldA, oldB); new bins to check + increment
                        // are everything in [newA, newB) outside [oldA, oldB).
                        double newRtStartPad = Math.Min(slot.RtStart, rtStartPad);
                        double newRtStopPad = Math.Max(slot.RtStop, rtStopPad);
                        var (oldA, oldB) = RtToBinRange(slot.RtStart, slot.RtStop);
                        var (extA, extB) = RtToBinRange(newRtStartPad, newRtStopPad);
                        int extLoFrom = extA, extLoTo = Math.Min(oldA, extB);
                        int extHiFrom = Math.Max(oldB, extA), extHiTo = extB;
                        if (extLoFrom < extLoTo && RangeMax(slotsPerBin, extLoFrom, extLoTo) + 1 > parameters.CycleBudget)
                            continue;
                        if (extHiFrom < extHiTo && RangeMax(slotsPerBin, extHiFrom, extHiTo) + 1 > parameters.CycleBudget)
                            continue;
                        if (OverlapsOtherSlotBand(sid, newMzMin - prmHalf, newMzMax + prmHalf,
                                                  newRtStartPad, newRtStopPad))
                            continue;
                        if (extLoFrom < extLoTo) RangeAddOne(slotsPerBin, extLoFrom, extLoTo);
                        if (extHiFrom < extHiTo) RangeAddOne(slotsPerBin, extHiFrom, extHiTo);
                        slot.RtStart = newRtStartPad;
                        slot.RtStop = newRtStopPad;
                        slot.CoStart = newCoStart;
                        slot.CoStop = newCoStop;
                        slot.MzMin = newMzMin;
                        slot.MzMax = newMzMax;
                        var merged = new double[slot.Fragments.Length + frags.Length];
                        Array.Copy(slot.Fragments, merged, slot.Fragments.Length);
                        Array.Copy(frags, 0, merged, slot.Fragments.Length, frags.Length);
                        Array.Sort(merged);
                        slot.Fragments = merged;
                        slot.MemberIndices.Add(candIndex);
                        RegisterSlotMzBins(sid, newMzMin, newMzMax);
                        TightenBounds(slot, candIndex);
                        return sid;
                    }
                }
            }

            // No existing slot accepted - open a new one if there's budget.
            var (newA, newB) = RtToBinRange(rtStartPad, rtStopPad);
            if (newA >= newB) return -1;
            if (RangeMax(slotsPerBin, newA, newB) + 1 > parameters.CycleBudget) return -1;
            // A fresh solo band that would already overlap a concurrent slot's
            // band cannot be placed disjointly; the caller falls through to
            // this protein's next candidate.
            if (OverlapsOtherSlotBand(-1, mz - prmHalf, mz + prmHalf, rtStartPad, rtStopPad)) return -1;
            RangeAddOne(slotsPerBin, newA, newB);
            int newSid = slots.Count;
            var newSlot = new Slot
            {
                Id = newSid,
                MzMin = mz,
                MzMax = mz,
                RtStart = rtStartPad,
                RtStop = rtStopPad,
                CoStart = rtStart,
                CoStop = rtStop,
                Fragments = frags.Length == 0 ? Array.Empty<double>() : (double[])frags.Clone(),
            };
            newSlot.MemberIndices.Add(candIndex);
            slots.Add(newSlot);
            RegisterSlotMzBins(newSid, mz, mz);
            TightenBounds(newSlot, candIndex);
            return newSid;
        }

        // Group candidates by protein, then optionally filter by target list.
        var groupToRanked = new Dictionary<string, List<int>>();
        for (int i = 0; i < candidates.Count; i++)
        {
            var c = candidates[i];
            if (!groupToRanked.TryGetValue(c.ProteinGroup, out var list))
            {
                list = new List<int>();
                groupToRanked[c.ProteinGroup] = list;
            }
            list.Add(i);
        }

        bool hasTargets = parameters.TargetProteins.Count > 0;

        // Match a candidate's protein group against the user's target set.
        // Skyline-document ingest carries the full Skyline protein name
        // (e.g. "sp|P55011|S12A2_HUMAN") as the group identifier, while
        // ProteinListParser extracts bare accessions ("P55011") from FASTA
        // headers. Without normalisation Exclusive mode silently drops
        // every group whose name carries the pipe-prefix, so we also try
        // the extracted accession on a miss.
        bool IsTargetGroup(string group)
        {
            if (parameters.TargetProteins.Contains(group)) return true;
            string accession = Ingest.ProteinListParser.ExtractAccession(group);
            return !ReferenceEquals(accession, group)
                && parameters.TargetProteins.Contains(accession);
        }

        if (hasTargets && parameters.TargetMode == TargetListMode.Exclusive)
        {
            var keep = new Dictionary<string, List<int>>();
            foreach (var (group, idxs) in groupToRanked)
            {
                if (IsTargetGroup(group))
                    keep[group] = idxs;
            }
            groupToRanked = keep;
        }

        // Rank two peptides within a group: unique-before-razor tier, then the
        // given metric (descending quantity / predicted score, ascending
        // q-value), then PrecursorId for a stable TOTAL order. The ranking
        // metric is a PARAMETER, not fixed to parameters.PeptideRanking, so the
        // tier loop can rank the confirmed and predicted tiers by different
        // metrics - see the ranking selection in that loop.
        //
        // The PrecursorId tiebreak matters: equal-scoring peptides would
        // otherwise keep whatever relative order the (unstable) sort left them
        // in, which depends on the order the candidate list arrived in - and
        // that order is NOT stable (CarafeLibraryReader materialises candidates
        // in Dictionary enumeration order, and .NET randomizes string hashing
        // per process). Ordering by the unique PrecursorId makes tie resolution
        // reproducible across app launches.
        int RankCompare(int x, int y, PeptidePriority ranking)
        {
            var cx = candidates[x];
            var cy = candidates[y];
            int tierX = cx.PeptideType == "unique" ? 0 : cx.PeptideType == "razor" ? 1 : 2;
            int tierY = cy.PeptideType == "unique" ? 0 : cy.PeptideType == "razor" ? 1 : 2;
            if (tierX != tierY) return tierX - tierY;
            int cmp = ranking switch
            {
                PeptidePriority.QValue => cx.QValue.CompareTo(cy.QValue),
                PeptidePriority.PredictedRanking =>
                    (cy.PredictedPriority ?? double.NegativeInfinity)
                        .CompareTo(cx.PredictedPriority ?? double.NegativeInfinity),
                _ => cy.PrecursorQuantity.CompareTo(cx.PrecursorQuantity),
            };
            if (cmp != 0) return cmp;
            return string.CompareOrdinal(cx.PrecursorId, cy.PrecursorId);
        }

        // Target-list rank is tier-independent (a group is on the user's target
        // list or it isn't). The rest of the group ordering - priority metric
        // and the smallest-group tiebreaker - is computed PER TIER inside the
        // scheduling loop below, over that tier's candidates only, so the
        // experimental tier orders proteins exactly as it would with no
        // predicted pool present.
        int TargetRank(string group)
        {
            if (!hasTargets) return 0;
            return IsTargetGroup(group) ? 0 : 1;
        }

        var allGroups = groupToRanked.Keys.ToList();

        var scheduledToSlot = new Dictionary<int, int>();
        var scheduledPerGroup = new Dictionary<string, int>(allGroups.Count);
        foreach (var g in allGroups) { scheduledPerGroup[g] = 0; }
        var coveredGroups = new HashSet<string>();

        // Shrink-with-re-absorption (WindowPacking.Disjoint only).
        //
        // Greedy slots only ever GROW, so a precursor that cannot merge and
        // cannot open a disjoint window has nowhere to go and is dropped. This
        // recovers a slice of those: if the blocking slot has members far
        // enough away to survive, evict the ones inside the new precursor's
        // quad width and let the NEW window absorb them as multiplexed members.
        // They sit within PrmIsolationWidthTh of it, so they are already in its
        // m/z neighbourhood; most co-elute and don't clash, and ride the same
        // chimeric scan. Nothing is dropped unless it cannot be re-absorbed AND
        // its protein would survive the loss.
        int TryShrinkAndPlace(int candIndex)
        {
            if (!disjointPacking) return -1;
            var prec = candidates[candIndex];
            if (prec.Top4Fragments.Length == 0) return -1;

            double mz = prec.PrecursorMz;
            double soloLo = mz - prmHalf, soloHi = mz + prmHalf;
            double padStart = prec.RtStart - padMin, padStop = prec.RtStop + padMin;
            double centerSpanBudget =
                Math.Max(0.0, parameters.IsolationWindowTh - parameters.PrmIsolationWidthTh);

            // Concurrent slots whose band overlaps the prospective solo band.
            var blockers = new List<int>();
            for (int b = (int)Math.Floor(soloLo) - 1; b <= (int)Math.Floor(soloHi) + 1; b++)
            {
                if (!slotByMzBin.TryGetValue(b, out var sids)) continue;
                foreach (int sid in sids)
                {
                    if (blockers.Contains(sid)) continue;
                    var s = slots[sid];
                    if (!(padStart < s.RtStop && s.RtStart < padStop)) continue;
                    if (soloLo < s.MzMax + prmHalf && s.MzMin - prmHalf < soloHi) blockers.Add(sid);
                }
            }
            if (blockers.Count == 0) return -1;

            // Locked experimental tier: a predicted candidate may never shrink,
            // evict, or re-geometry a slot that holds an experimental member.
            // Shrink only ever rearranges predicted-only slots. During the
            // experimental tier the placing candidate is itself experimental, so
            // this is a no-op and the published shrink behaviour (and the golden
            // fixtures) are unchanged.
            if (prec.Source == CandidateSource.Predicted)
            {
                foreach (int sid in blockers)
                    if (slots[sid].MemberIndices.Any(
                            m => candidates[m].Source == CandidateSource.Experimental))
                        return -1;
            }

            // Dry-run the shrink. Every blocker must keep at least one member a
            // full quad width away, otherwise "shrinking" it means emptying it.
            var evicted = new List<int>();
            var shrunk = new Dictionary<int, (List<int> Kept, double MzMin, double MzMax,
                                              double Co0, double Co1, double P0, double P1)>();
            foreach (int sid in blockers)
            {
                var kept = new List<int>();
                foreach (int m in slots[sid].MemberIndices)
                {
                    if (Math.Abs(candidates[m].PrecursorMz - mz) >= parameters.PrmIsolationWidthTh)
                        kept.Add(m);
                    else
                        evicted.Add(m);
                }
                if (kept.Count == 0) return -1;

                double kMin = double.PositiveInfinity, kMax = double.NegativeInfinity;
                double kCo0 = double.NegativeInfinity, kCo1 = double.PositiveInfinity;
                double kP0 = double.PositiveInfinity, kP1 = double.NegativeInfinity;
                foreach (int m in kept)
                {
                    var mc = candidates[m];
                    kMin = Math.Min(kMin, mc.PrecursorMz); kMax = Math.Max(kMax, mc.PrecursorMz);
                    kCo0 = Math.Max(kCo0, mc.RtStart); kCo1 = Math.Min(kCo1, mc.RtStop);
                    kP0 = Math.Min(kP0, mc.RtStart - padMin); kP1 = Math.Max(kP1, mc.RtStop + padMin);
                }
                shrunk[sid] = (kept, kMin, kMax, kCo0, kCo1, kP0, kP1);
            }
            if (evicted.Count == 0) return -1;

            // Re-absorb what we can into the new window.
            var members = new List<int> { candIndex };
            var homeless = new List<int>();
            double winLo = mz, winHi = mz, coStart = prec.RtStart, coStop = prec.RtStop;
            double winPadStart = padStart, winPadStop = padStop;
            foreach (int e in evicted)
            {
                var ec = candidates[e];
                double nLo = Math.Min(winLo, ec.PrecursorMz), nHi = Math.Max(winHi, ec.PrecursorMz);
                double nCo0 = Math.Max(coStart, ec.RtStart), nCo1 = Math.Min(coStop, ec.RtStop);
                bool fits =
                    (nHi - nLo) <= centerSpanBudget
                    && nCo0 < nCo1
                    && !members.Any(w => FragmentClash.AnyWithin(
                           ec.Top4Fragments, candidates[w].Top4Fragments, parameters.FragmentTolDa))
                    && (parameters.ChargeHandling != ChargeHandling.SameChargePerSlot
                        || ec.PrecursorCharge == prec.PrecursorCharge);
                if (fits)
                {
                    members.Add(e);
                    winLo = nLo; winHi = nHi; coStart = nCo0; coStop = nCo1;
                    winPadStart = Math.Min(winPadStart, ec.RtStart - padMin);
                    winPadStop = Math.Max(winPadStop, ec.RtStop + padMin);
                }
                else homeless.Add(e);
            }

            // Never un-cover a protein to make room.
            foreach (int h in homeless)
                if (scheduledPerGroup.TryGetValue(candidates[h].ProteinGroup, out int n) && n <= 1)
                    return -1;

            // The new window is WIDER than the solo band that found the
            // blockers, so re-verify the final geometry against every slot -
            // blockers at their post-shrink size, everything else as-is.
            double winBandLo = winLo - prmHalf, winBandHi = winHi + prmHalf;
            for (int sid = 0; sid < slots.Count; sid++)
            {
                double sMin, sMax, sRt0, sRt1;
                if (shrunk.TryGetValue(sid, out var g))
                {
                    sMin = g.MzMin; sMax = g.MzMax; sRt0 = g.P0; sRt1 = g.P1;
                }
                else
                {
                    var s = slots[sid];
                    if (s.MemberIndices.Count == 0) continue;
                    sMin = s.MzMin; sMax = s.MzMax; sRt0 = s.RtStart; sRt1 = s.RtStop;
                }
                if (!(winPadStart < sRt1 && sRt0 < winPadStop)) continue;
                if (winBandLo < sMax + prmHalf && sMin - prmHalf < winBandHi) return -1;
            }

            // Budget for the new window, checked before mutating anything.
            var (winA, winB) = RtToBinRange(winPadStart, winPadStop);
            if (winA >= winB) return -1;
            if (RangeMax(slotsPerBin, winA, winB) + 1 > parameters.CycleBudget) return -1;

            // ---- commit ----
            foreach (var (sid, g) in shrunk)
            {
                var s = slots[sid];
                var (oldA, oldB) = RtToBinRange(s.RtStart, s.RtStop);

                s.MemberIndices.Clear();
                var frags = new List<double>();
                foreach (int m in g.Kept)
                {
                    s.MemberIndices.Add(m);
                    frags.AddRange(candidates[m].Top4Fragments);
                }
                frags.Sort();
                s.MzMin = g.MzMin; s.MzMax = g.MzMax;
                s.CoStart = g.Co0; s.CoStop = g.Co1;
                s.RtStart = g.P0; s.RtStop = g.P1;
                s.Fragments = frags.ToArray();
                // Growth bounds were derived from members that may now be gone.
                // Leaving them is safe: they only ever tighten, so a stale bound
                // is conservative, never permissive.

                // Release the occupancy the slot no longer spans.
                var (newA2, newB2) = RtToBinRange(s.RtStart, s.RtStop);
                for (int i = oldA; i < oldB; i++)
                    if (i < newA2 || i >= newB2) slotsPerBin[i]--;
            }

            RangeAddOne(slotsPerBin, winA, winB);
            int sidNew = slots.Count;
            var win = new Slot
            {
                Id = sidNew,
                MzMin = winLo,
                MzMax = winHi,
                RtStart = winPadStart,
                RtStop = winPadStop,
                CoStart = coStart,
                CoStop = coStop,
            };
            var wf = new List<double>();
            foreach (int m in members) { win.MemberIndices.Add(m); wf.AddRange(candidates[m].Top4Fragments); }
            wf.Sort();
            win.Fragments = wf.ToArray();
            slots.Add(win);
            RegisterSlotMzBins(sidNew, winLo, winHi);
            foreach (int m in members) TightenBounds(win, m);

            // Re-absorbed members changed slots but stay scheduled; homeless
            // ones lose their place (their protein keeps another peptide).
            foreach (int m in members) if (m != candIndex) scheduledToSlot[m] = sidNew;
            foreach (int h in homeless)
            {
                scheduledToSlot.Remove(h);
                string g = candidates[h].ProteinGroup;
                if (scheduledPerGroup.TryGetValue(g, out int n)) scheduledPerGroup[g] = n - 1;
            }
            return sidNew;
        }

        // Place a candidate: the normal merge-or-open, then the shrink fallback.
        int Place(int candIndex)
        {
            int sid = TrySchedule(candIndex);
            return sid >= 0 ? sid : TryShrinkAndPlace(candIndex);
        }
        // Effective load-up cap by objective:
        //   Balanced         - user's MaxPeptidesPerProtein.
        //   MaximizeProteins - user's MaxPeptidesPerProtein. The cover
        //                      pass is single-pass, so once coverage
        //                      converges there's no future cover pass
        //                      to save budget for. Capping load-up
        //                      below Max would leave cycle budget on
        //                      the table after every protein had had
        //                      its first peptide; load-up should fill
        //                      to Max instead. The difference vs
        //                      Balanced is the cover-pass strategy
        //                      (prefer-joinable), not the load-up cap.
        //   MaximizePeptides - no per-group cap; round-robin order in
        //                      the load-up loop preserves fairness.
        int rawMax = Math.Max(1, parameters.MaxPeptidesPerProtein);
        int effectiveMaxPerGroup = parameters.Objective switch
        {
            CoverageObjective.MaximizePeptides => int.MaxValue,
            _ => rawMax,
        };

        // Partition each group's ranked queue by source tier, preserving the
        // within-group sort. Predicted candidates land in a separate queue so
        // they only compete for budget after the experimental tier is done.
        // Per-tier queues, each sorted by that tier's own ranking metric. The
        // confirmed tier is ranked by real evidence (quantity / q-value) and the
        // predicted tier by its predicted score - see the ranking selection in
        // the tier loop.
        Dictionary<string, List<int>> QueuesForTier(CandidateSource src, PeptidePriority ranking) =>
            allGroups.ToDictionary(
                g => g,
                g =>
                {
                    var list = groupToRanked[g].Where(i => candidates[i].Source == src).ToList();
                    list.Sort((x, y) => RankCompare(x, y, ranking));
                    return list;
                });

        bool useLookAhead = parameters.Objective != CoverageObjective.Balanced;

        // True when some Experimental-tagged candidate carries a predicted score.
        // This is the Carafe-ONLY pool (predicted peptides tagged Experimental by
        // default, WITH a detectability / Bromo score): there PredictedRanking is
        // meaningful for the experimental tier. In the DIA-NN + Carafe backfill
        // pool the confirmed peptides have no score, so the experimental tier must
        // not rank by PredictedRanking. Distinguishes the two.
        bool experimentalHasPredictedScore = false;
        for (int i = 0; i < candidates.Count; i++)
        {
            if (candidates[i].Source == CandidateSource.Experimental
                && candidates[i].PredictedPriority is not null)
            {
                experimentalHasPredictedScore = true;
                break;
            }
        }

        // Two-tier scheduling. Schedule every Experimental (measured) candidate
        // first via the full cover + load-up, then let Predicted (Carafe
        // backfill) candidates reuse whatever cycle budget is left. The shared
        // occupancy (slotsPerBin / slots) and per-group counts persist across
        // tiers, so predicted peptides respect the same budget and the same
        // per-protein cap as the confirmed ones. A single-tier all-experimental
        // pool (every existing ingest path and test) leaves the predicted
        // queues empty, so the second iteration is a no-op and the schedule is
        // identical to before - this is what preserves the golden fixtures.
        foreach (var tier in new[] { CandidateSource.Experimental, CandidateSource.Predicted })
        {
            // Tier-aware peptide ranking. The two tiers rank by DIFFERENT metrics:
            //
            //  - Predicted tier: always its predicted score (Carafe detectability
            //    / Bromo). This is independent of the confirmed tier's metric, so
            //    enabling backfill never has to change - and thus never disturbs -
            //    the confirmed peptides' ranking.
            //  - Confirmed tier: the user's chosen evidence metric
            //    (PrecursorQuantity or QValue), so the confirmed schedule is
            //    identical to a backfill-off run. If that metric is somehow
            //    PredictedRanking while the confirmed peptides carry no predicted
            //    score (real DIA-NN targets), fall back to PrecursorQuantity so
            //    they aren't collapsed to a PrecursorId tie. A Carafe-only pool
            //    tags its predicted peptides Experimental but DOES give them
            //    scores, so it keeps PredictedRanking (experimentalHasPredictedScore).
            var tierRanking =
                tier == CandidateSource.Predicted
                    ? PeptidePriority.PredictedRanking
                : parameters.PeptideRanking == PeptidePriority.PredictedRanking
                        && !experimentalHasPredictedScore
                    ? PeptidePriority.PrecursorQuantity
                    : parameters.PeptideRanking;
            var queues = QueuesForTier(tier, tierRanking);
            var cursors = allGroups.ToDictionary(g => g, _ => 0);

            // Rebuild the disjoint-packing clash-barrier pool for this tier.
            // Experimental tier: only experimental precursors can bound an
            // experimental slot's growth, so the bounds (and the schedule) don't
            // depend on the predicted pool. Predicted tier: every precursor is a
            // real interferer for a predicted slot, so the pool is the full set.
            if (disjointPacking)
            {
                candByMzBin.Clear();
                for (int i = 0; i < candidates.Count; i++)
                {
                    if (tier == CandidateSource.Experimental
                        && candidates[i].Source != CandidateSource.Experimental) continue;
                    int bin = (int)Math.Floor(candidates[i].PrecursorMz);
                    if (!candByMzBin.TryGetValue(bin, out var list))
                    {
                        list = new List<int>();
                        candByMzBin[bin] = list;
                    }
                    list.Add(i);
                }
            }

            // Per-tier group ordering. Priority metric (summed intensity / best
            // q-value) and the smallest-group tiebreaker are computed over THIS
            // tier's candidates only, so the experimental tier visits proteins in
            // exactly the order it would with no predicted pool present. The
            // predicted tier ranks over its own candidates for the same reason it
            // schedules last: it only ever fills leftover budget.
            var tierSummed = new Dictionary<string, double>(allGroups.Count);
            var tierBestQ = new Dictionary<string, double>(allGroups.Count);
            foreach (var g in allGroups)
            {
                double summed = 0.0;
                double bestQ = double.PositiveInfinity;
                foreach (var i in queues[g])
                {
                    summed += candidates[i].PrecursorQuantity;
                    double pq = candidates[i].ProteinQValue;
                    if (!double.IsNaN(pq) && pq < bestQ) bestQ = pq;
                }
                tierSummed[g] = summed;
                tierBestQ[g] = bestQ;
            }

            int PriorityRank(string group) => parameters.ProteinRanking switch
            {
                ProteinPriority.SummedIntensity =>
                    -(int)Math.Min(tierSummed[group] / 1e3, int.MaxValue),
                ProteinPriority.ProteinQValue =>
                    (int)Math.Min(tierBestQ[group] * 1e6, int.MaxValue),
                // Without a positional input we treat it as alphabetical; the UI
                // passes targets as an ordered list and that order is preserved
                // by the in-target / out-of-target split.
                ProteinPriority.ProvidedListOrder => 0,
                _ => 0,
            };

            var groupOrder = allGroups.ToList();
            groupOrder.Sort((g1, g2) =>
            {
                int t = TargetRank(g1).CompareTo(TargetRank(g2));
                if (t != 0) return t;
                int p = PriorityRank(g1).CompareTo(PriorityRank(g2));
                if (p != 0) return p;
                int c = queues[g1].Count.CompareTo(queues[g2].Count);
                if (c != 0) return c;
                return string.CompareOrdinal(g1, g2);
            });

        // Pass 1: cover one peptide per group. Cover-pass strategy is
        // selected by Objective:
        //
        // - Balanced (published webinar): reactive cover. Walk the
        //   static intensity-sorted queue and TrySchedule; on a per-bin
        //   budget hit, fall back to the next peptide in score order.
        //   Slide 9 exactly.
        //
        // - MaximizeProteins: look-ahead cover with prefer-joinable.
        //   First walk the queue in score order and accept any peptide
        //   that joins an existing slot for free; if none joinable,
        //   fall through to the look-ahead pick (least-saturated RT
        //   bin) and finally to a static-order fallback. Conserves
        //   budget for first-peptide coverage of more proteins.
        //
        // - MaximizePeptides: look-ahead cover only (no joinable-first
        //   preference); same load-up shape as MaximizeProteins's
        //   look-ahead path, but with no cap on load-up.
        //
        // The Objective also decides the cursor advancement: reactive
        // cover advances the cursor past peptides it tried so they
        // don't get retried in the load-up; look-ahead cover holds the
        // cursor at 0 because it considered-but-didn't-try most of the
        // queue (load-up walks from 0 and skips already-scheduled
        // candidates to recover them).
        foreach (var g in groupOrder)
        {
            cancellationToken.ThrowIfCancellationRequested();
            // Skip groups an earlier tier already covered: a predicted peptide
            // never re-covers a protein an experimental peptide already
            // covered; its depth is added by the load-up loop instead.
            if (coveredGroups.Contains(g)) continue;
            var queue = queues[g];
            int cursor = cursors[g];
            int chosenSid = -1;
            int chosenCandIdx = -1;

            // MaximizeProteins: prefer a peptide that joins an existing
            // slot for free over one that opens a new slot. Walk the queue
            // in score order, accept the first joinable hit. If no peptide
            // is joinable, fall through to the regular cover branch so the
            // protein still gets covered by opening a slot (matching the
            // user's intent that this objective should not silently drop
            // hard-to-place proteins).
            if (parameters.Objective == CoverageObjective.MaximizeProteins)
            {
                for (int q = cursor; q < queue.Count; q++)
                {
                    int candIdx = queue[q];
                    if (!CanJoinExistingSlot(candIdx)) continue;
                    int sid = Place(candIdx);
                    if (sid >= 0)
                    {
                        chosenSid = sid;
                        chosenCandIdx = candIdx;
                        break;
                    }
                }
            }

            if (chosenSid < 0 && useLookAhead)
            {
                int bestCandIdx = -1;
                int bestMaxLoad = int.MaxValue;
                int bestTieRank = int.MaxValue;
                for (int q = cursor; q < queue.Count; q++)
                {
                    int candIdx = queue[q];
                    var c = candidates[candIdx];
                    double padStart = c.RtStart - padMin;
                    double padStop = c.RtStop + padMin;
                    var (a, b) = RtToBinRange(padStart, padStop);
                    if (a >= b) continue;
                    int load = RangeMax(slotsPerBin, a, b);
                    if (load >= parameters.CycleBudget) continue;
                    if (load < bestMaxLoad
                        || (load == bestMaxLoad && q < bestTieRank))
                    {
                        bestMaxLoad = load;
                        bestCandIdx = candIdx;
                        bestTieRank = q;
                    }
                }
                if (bestCandIdx >= 0)
                {
                    int sid = Place(bestCandIdx);
                    if (sid >= 0)
                    {
                        chosenSid = sid;
                        chosenCandIdx = bestCandIdx;
                    }
                }
                // Fallback if the RT-aware pick didn't fit (fragment
                // clash, co-elution, charge): walk the queue in static
                // order until something works.
                if (chosenSid < 0)
                {
                    for (int q = cursor; q < queue.Count; q++)
                    {
                        int candIdx = queue[q];
                        if (candIdx == bestCandIdx) continue;
                        int sid = Place(candIdx);
                        if (sid >= 0)
                        {
                            chosenSid = sid;
                            chosenCandIdx = candIdx;
                            break;
                        }
                    }
                }
                // Cursor stays at 0. Load-up will re-walk and skip
                // already-scheduled candidates.
            }
            else if (chosenSid < 0)
            {
                while (cursor < queue.Count)
                {
                    int candIdx = queue[cursor];
                    int sid = Place(candIdx);
                    cursor++;
                    if (sid >= 0)
                    {
                        chosenSid = sid;
                        chosenCandIdx = candIdx;
                        break;
                    }
                }
                cursors[g] = cursor;
            }

            if (chosenSid >= 0)
            {
                scheduledToSlot[chosenCandIdx] = chosenSid;
                coveredGroups.Add(g);
                scheduledPerGroup[g]++;
            }
        }

        // Pass 2+: load-up loop. Breadth-first round-robin across
        // proteins: each outer iteration is one "lap" through every
        // covered protein, and each protein adds AT MOST ONE peptide
        // per lap. So after K laps every protein has been given K
        // chances at adding its next peptide. A protein drops out of
        // the loop when it either hits its cap or its candidate queue
        // runs dry. The outer while exits when a full lap adds
        // nothing.
        //
        // Per-Objective cap (effectiveMaxPerGroup, computed above):
        //   Balanced         -> cap = MaxPeptidesPerProtein (default 5).
        //   MaximizeProteins -> cap = MaxPeptidesPerProtein (default 5).
        //   MaximizePeptides -> cap = int.MaxValue.
        {
            while (true)
            {
                cancellationToken.ThrowIfCancellationRequested();
                bool changed = false;
                foreach (var g in groupOrder)
                {
                    if (!coveredGroups.Contains(g)) continue;
                    if (scheduledPerGroup[g] >= effectiveMaxPerGroup) continue;
                    var queue = queues[g];
                    int cursor = cursors[g];
                    while (cursor < queue.Count)
                    {
                        int candIdx = queue[cursor];
                        cursor++;
                        // Skip peptides already taken by cover pass.
                        // RT-aware cover may have chosen any position in
                        // the queue without advancing the cursor.
                        if (scheduledToSlot.ContainsKey(candIdx)) continue;
                        int sid = Place(candIdx);
                        if (sid >= 0)
                        {
                            scheduledToSlot[candIdx] = sid;
                            scheduledPerGroup[g]++;
                            changed = true;
                            break;
                        }
                    }
                    cursors[g] = cursor;
                }
                if (!changed) break;
            }
        }
        } // end tier loop (Experimental, then Predicted)

        // Min-peptides-per-protein enforcement: drop any group that
        // didn't reach the minimum and free the slot capacity its peptides
        // were holding. This runs once (not iterated) - in practice this
        // produces the user's expected "min=2 means every kept group has
        // >=2 peptides" semantics.
        if (parameters.MinPeptidesPerProtein > 1)
        {
            var groupCounts = new Dictionary<string, int>();
            foreach (var kv in scheduledToSlot)
            {
                string g = candidates[kv.Key].ProteinGroup;
                groupCounts[g] = (groupCounts.TryGetValue(g, out int n) ? n : 0) + 1;
            }
            var failedGroups = new HashSet<string>();
            foreach (var kv in groupCounts)
                if (kv.Value < parameters.MinPeptidesPerProtein) failedGroups.Add(kv.Key);

            if (failedGroups.Count > 0)
            {
                // Drop scheduled entries belonging to failed groups.
                var keptScheduled = new Dictionary<int, int>();
                foreach (var kv in scheduledToSlot)
                    if (!failedGroups.Contains(candidates[kv.Key].ProteinGroup))
                        keptScheduled[kv.Key] = kv.Value;
                scheduledToSlot = keptScheduled;
                foreach (var g in failedGroups) coveredGroups.Remove(g);

                // Group the kept entries by their (old) slot id so we can
                // rebuild slots with only the kept members and recompute
                // each slot's bounding box, co-elution window and fragment
                // set from scratch.
                var keptBySlot = new Dictionary<int, List<int>>();
                foreach (var kv in scheduledToSlot)
                {
                    if (!keptBySlot.TryGetValue(kv.Value, out var list))
                    {
                        list = new List<int>();
                        keptBySlot[kv.Value] = list;
                    }
                    list.Add(kv.Key);
                }

                var rebuiltSlots = new List<Slot>();
                var oldToNewSid = new Dictionary<int, int>();
                foreach (var oldSlot in slots)
                {
                    if (!keptBySlot.TryGetValue(oldSlot.Id, out var members)) continue;
                    int newSid = rebuiltSlots.Count;
                    oldToNewSid[oldSlot.Id] = newSid;
                    var fresh = new Slot { Id = newSid };
                    fresh.MemberIndices.AddRange(members);
                    double mzMin = double.PositiveInfinity, mzMax = double.NegativeInfinity;
                    double coStart = double.NegativeInfinity, coStop = double.PositiveInfinity;
                    double padStart = double.PositiveInfinity, padStop = double.NegativeInfinity;
                    var mergedFrags = new List<double>();
                    foreach (var idx in members)
                    {
                        var c = candidates[idx];
                        if (c.PrecursorMz < mzMin) mzMin = c.PrecursorMz;
                        if (c.PrecursorMz > mzMax) mzMax = c.PrecursorMz;
                        if (c.RtStart > coStart) coStart = c.RtStart;
                        if (c.RtStop < coStop) coStop = c.RtStop;
                        double ps = c.RtStart - padMin, pe = c.RtStop + padMin;
                        if (ps < padStart) padStart = ps;
                        if (pe > padStop) padStop = pe;
                        mergedFrags.AddRange(c.Top4Fragments);
                    }
                    fresh.MzMin = mzMin;
                    fresh.MzMax = mzMax;
                    fresh.CoStart = coStart;
                    fresh.CoStop = coStop;
                    fresh.RtStart = padStart;
                    fresh.RtStop = padStop;
                    mergedFrags.Sort();
                    fresh.Fragments = mergedFrags.ToArray();
                    rebuiltSlots.Add(fresh);
                }

                // Remap scheduled slot ids to the rebuilt slot list.
                foreach (var key in scheduledToSlot.Keys.ToList())
                    scheduledToSlot[key] = oldToNewSid[scheduledToSlot[key]];
                slots = rebuiltSlots;

                // Recompute the concurrent-slot occupancy curve from the
                // rebuilt slots so peak load reflects the surviving
                // schedule (not the budget that was spent and then freed).
                Array.Clear(slotsPerBin);
                foreach (var slot in slots)
                {
                    var (a, b) = RtToBinRange(slot.RtStart, slot.RtStop);
                    for (int i = a; i < b; i++) slotsPerBin[i]++;
                }
            }
        }

        var scheduledIndices = scheduledToSlot.Keys.ToArray();
        var scheduledSlotIds = new int[scheduledIndices.Length];
        for (int i = 0; i < scheduledIndices.Length; i++)
            scheduledSlotIds[i] = scheduledToSlot[scheduledIndices[i]];

        var rtGrid = new double[nBins];
        for (int i = 0; i < nBins; i++)
            rtGrid[i] = rtLo + (i + 0.5) * parameters.RtBinMin;

        int proteinsCovered = coveredGroups.Count;

        return new ScheduleResult
        {
            ScheduledIndices = scheduledIndices,
            ScheduledSlotIds = scheduledSlotIds,
            Slots = slots.ToArray(),
            RtGrid = rtGrid,
            SlotCountCurve = slotsPerBin,
            ProteinGroupsCovered = proteinsCovered,
        };
    }
}
