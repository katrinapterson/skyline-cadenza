using SkylineCadenza.Core.Scheduling;
using Xunit;

namespace SkylineCadenza.Tests;

/// <summary>
/// Covers the window-scoping invariants described in
/// <c>docs/scheduler-window-scoping.md</c>. A slot's fragment-clash check is
/// scoped to its members while the isolation window physically co-isolates
/// everything in its m/z band, so without per-slot growth bounds expanding a
/// window could sweep it across a precursor that was explicitly rejected for
/// clashing (the "illusory exclusion"). Slot-aware packing - now the only
/// behaviour - prevents that with growth bounds and a neighbour-band probe.
/// </summary>
public class SchedulerWindowScopingTests
{
    private const double Prm = 0.7;
    private const double Tol = 0.05;

    private static Candidate Make(
        string id, double mz, double rtStart, double rtStop, string group,
        double quantity, double[] fragments)
    {
        var ions = new FragmentIon[fragments.Length];
        for (int i = 0; i < fragments.Length; i++)
            ions[i] = new FragmentIon(fragments[i], fragments.Length - i, 1);
        return new Candidate
        {
            PrecursorId = id,
            StrippedSequence = id,
            ModifiedSequence = id,
            PrecursorCharge = 2,
            PrecursorMz = mz,
            RtStart = rtStart,
            RtStop = rtStop,
            RtApex = (rtStart + rtStop) / 2.0,
            PrecursorQuantity = quantity,
            QValue = 0.001,
            ProteinQValue = 0.001,
            Proteotypic = 1,
            ProteinGroup = group,
            PeptideType = "unique",
            Fragments = ions,
            Top4Fragments = Candidate.DeriveTopMz(ions, 4),
            Run = "test",
        };
    }

    // A and P share no fragment (may multiplex); Q sits between them in m/z and
    // clashes with EACH. The A-P center span is 1.5 Th, inside the 2.3 Th
    // budget, so an unbounded merge would produce a band [599.65, 601.85] that
    // swallows Q (600.75) - the precursor it just rejected for clashing with A.
    // A and Q solo (0.75 Th apart, 0.7 Th wide) do NOT co-isolate; only the
    // expansion does, which is exactly what the growth bounds prevent.
    private static Candidate[] AqpCandidates() => new[]
    {
        Make("A", 600.00, 5.0, 5.4, "G1", 2e6, new[] { 1000.0, 1100.0, 1200.0, 1300.0 }),
        Make("Q", 600.75, 5.0, 5.4, "G2", 1e6, new[] { 1000.0, 1400.0, 1500.0, 1600.0 }),
        Make("P", 601.50, 5.0, 5.4, "G1", 1e6, new[] { 1400.0, 1700.0, 1800.0, 1900.0 }),
    };

    private static SchedulingParameters Params() => new()
    {
        CycleBudget = 100,
        IsolationWindowTh = 3.0,
        PrmIsolationWidthTh = Prm,
        FragmentTolDa = Tol,
        Objective = CoverageObjective.MaximizeProteins,
    };

    /// <summary>
    /// Counts slots whose isolation band co-isolates a co-eluting non-member
    /// that clashes with one of the slot's members.
    /// </summary>
    private static int CountUncheckedCoIsolations(Candidate[] cands, ScheduleResult result)
    {
        int violations = 0;
        foreach (var slot in result.Slots)
        {
            double bandLo = slot.MzMin - Prm / 2.0;
            double bandHi = slot.MzMax + Prm / 2.0;
            for (int ci = 0; ci < cands.Length; ci++)
            {
                if (slot.MemberIndices.Contains(ci)) continue;
                var c = cands[ci];
                if (c.PrecursorMz < bandLo || c.PrecursorMz > bandHi) continue;
                if (Math.Max(slot.CoStart, c.RtStart) >= Math.Min(slot.CoStop, c.RtStop)) continue;
                if (slot.MemberIndices.Any(mi =>
                        FragmentClash.AnyWithin(c.Top4Fragments, cands[mi].Top4Fragments, Tol)))
                    violations++;
            }
        }
        return violations;
    }

    [Fact]
    public void Scheduler_SlotAwarePacking_WindowNeverExpandsOverClashingNonMember()
    {
        var cands = AqpCandidates();
        var result = Scheduler.Run(cands, Params());

        Assert.Equal(0, CountUncheckedCoIsolations(cands, result));

        // A's growth ceiling must stop short of Q, so P cannot join A's slot.
        int aIdx = Array.FindIndex(cands, c => c.PrecursorId == "A");
        int pIdx = Array.FindIndex(cands, c => c.PrecursorId == "P");
        var aSlot = result.Slots.First(s => s.MemberIndices.Contains(aIdx));
        Assert.DoesNotContain(pIdx, aSlot.MemberIndices);
    }

    [Fact]
    public void Scheduler_PrmMode_SlotsStaySoloAndZeroSpan()
    {
        // A PRM slot is one precursor at a fixed quadrupole width, so the
        // slot-aware packing guards must not engage: no growth bounds, no
        // neighbour-band rejections, every slot solo and zero-span.
        var cands = AqpCandidates();
        var prm = Params() with { Mode = AcquisitionMode.Prm };

        var result = Scheduler.Run(cands, prm);

        Assert.Equal(3, result.ScheduledIndices.Length);
        Assert.All(result.Slots, s =>
        {
            Assert.Single(s.MemberIndices);
            Assert.Equal(s.MzMin, s.MzMax, 9);
        });
    }

    [Fact]
    public void Scheduler_SlotAwarePacking_SurvivingCoIsolationsAreAllUnseparable()
    {
        // Disjoint cannot eliminate every co-isolation: two co-eluting
        // precursors closer than PrmIsolationWidthTh cannot be given disjoint
        // bands by ANY window scheme, so those survive and can only be resolved
        // by dropping one. What must NOT survive is a separable one - an
        // interferer far enough away that a tighter band would have excluded it.
        //
        // s0/s1 clash 0.2 Th apart -> unseparable, may survive.
        //
        // s2/x clash 1.0 Th apart -> separable. A lone clashing pair is never
        // co-isolated by solo bands (they just take separate slots), so
        // the discriminating case needs a THIRD precursor that widens a band
        // across the interferer: s4 doesn't clash s2, so an unbounded merge
        // would join s2+s4
        // (span 2.0 Th, inside the 2.3 budget) into a band [601.65, 604.35] that
        // swallows x at 603.0. Disjoint must bound s2's growth short of x.
        var cands = new[]
        {
            Make("s0", 600.0, 5.0, 5.4, "G0", 6e6, new[] { 1000.0, 1100.0, 1200.0, 1300.0 }),
            Make("s1", 600.2, 5.0, 5.4, "G1", 5e6, new[] { 1000.0, 2100.0, 2200.0, 2300.0 }),
            Make("s2", 602.0, 5.0, 5.4, "G2", 4e6, new[] { 1000.0, 2000.0, 3100.0, 3200.0 }),
            Make("x",  603.0, 5.0, 5.4, "G3", 3e6, new[] { 2000.0, 4100.0, 4200.0, 4300.0 }),
            Make("s4", 604.0, 5.0, 5.4, "G4", 2e6, new[] { 5000.0, 5100.0, 5200.0, 5300.0 }),
        };

        var result = Scheduler.Run(cands, Params());

        foreach (var slot in result.Slots)
        {
            double bandLo = slot.MzMin - Prm / 2.0;
            double bandHi = slot.MzMax + Prm / 2.0;
            for (int ci = 0; ci < cands.Length; ci++)
            {
                if (slot.MemberIndices.Contains(ci)) continue;
                var c = cands[ci];
                if (c.PrecursorMz < bandLo || c.PrecursorMz > bandHi) continue;
                if (Math.Max(slot.CoStart, c.RtStart) >= Math.Min(slot.CoStop, c.RtStop)) continue;

                var clashing = slot.MemberIndices
                    .Where(mi => FragmentClash.AnyWithin(c.Top4Fragments, cands[mi].Top4Fragments, Tol))
                    .ToList();
                if (clashing.Count == 0) continue;

                double gap = clashing.Min(mi => Math.Abs(c.PrecursorMz - cands[mi].PrecursorMz));
                Assert.True(gap < Prm,
                    $"'{c.PrecursorId}' is co-isolated by a slot it clashes with from {gap:F2} Th away " +
                    $"(>= {Prm} Th), so a tighter band could have excluded it. Disjoint must only " +
                    "leave sub-quad-width co-isolations, which no window scheme can separate.");
            }
        }
    }

    /// <summary>
    /// A wide slot blocks a later precursor, but has a member far enough away
    /// to survive shrinking. The blocked precursor should take the space, and
    /// the evicted member should be RE-ABSORBED into the new window rather than
    /// dropped - it sits inside the new window's m/z neighbourhood anyway.
    ///
    ///   A 600.0 (G1) opens a slot.
    ///   B 602.0 (G2) joins it (span 2.0, inside the 2.3 budget) -> band [599.65, 602.35].
    ///   P 602.6 (G3) is blocked: its solo band [602.25, 602.95] overlaps that.
    ///
    /// Shrinking evicts B (0.6 Th from P, inside the quad width); A survives at
    /// 2.6 Th away. B then joins P's window, so all three proteins are covered
    /// and nothing is lost.
    /// </summary>
    [Fact]
    public void Scheduler_SlotAwarePacking_ShrinksBlockingSlotAndReabsorbsEvictedMember()
    {
        var cands = new[]
        {
            Make("A", 600.0, 5.0, 5.4, "G1", 3e6, new[] { 1000.0, 1100.0, 1200.0, 1300.0 }),
            Make("B", 602.0, 5.0, 5.4, "G2", 2e6, new[] { 2000.0, 2100.0, 2200.0, 2300.0 }),
            Make("P", 602.6, 5.0, 5.4, "G3", 1e6, new[] { 3000.0, 3100.0, 3200.0, 3300.0 }),
        };

        var result = Scheduler.Run(cands, Params());

        // All three proteins covered - P got in without costing B its place.
        Assert.Equal(3, result.ProteinGroupsCovered);
        Assert.Equal(3, result.ScheduledIndices.Length);

        int a = Array.FindIndex(cands, c => c.PrecursorId == "A");
        int b = Array.FindIndex(cands, c => c.PrecursorId == "B");
        int pi = Array.FindIndex(cands, c => c.PrecursorId == "P");

        // B was re-absorbed: it now shares P's window, not A's.
        var pSlot = result.Slots.First(s => s.MemberIndices.Contains(pi));
        Assert.Contains(b, pSlot.MemberIndices);
        var aSlot = result.Slots.First(s => s.MemberIndices.Contains(a));
        Assert.DoesNotContain(b, aSlot.MemberIndices);

        // A's slot shrank back to a solo band, and nothing overlaps.
        Assert.Single(aSlot.MemberIndices);
        AssertNoConcurrentBandOverlap(result);
    }

    /// <summary>
    /// The mirror case: the blocking slot has no member a full quad width away,
    /// so shrinking it would mean emptying it. The scheduler must refuse and
    /// leave the incumbent alone rather than trading one precursor for another.
    /// </summary>
    [Fact]
    public void Scheduler_SlotAwarePacking_RefusesToShrinkWhenBlockerCannotSurvive()
    {
        var cands = new[]
        {
            Make("A", 600.0, 5.0, 5.4, "G1", 3e6, new[] { 1000.0, 1100.0, 1200.0, 1300.0 }),
            // Clashes with A (1000) so it cannot merge, and sits 0.4 Th away so
            // A would have to evict its only member to make room.
            Make("P", 600.4, 5.0, 5.4, "G2", 1e6, new[] { 1000.0, 3100.0, 3200.0, 3300.0 }),
        };

        var result = Scheduler.Run(cands, Params());

        int a = Array.FindIndex(cands, c => c.PrecursorId == "A");
        int pi = Array.FindIndex(cands, c => c.PrecursorId == "P");
        Assert.Contains(a, result.ScheduledIndices);
        Assert.DoesNotContain(pi, result.ScheduledIndices);
        Assert.Equal(1, result.ProteinGroupsCovered);
        AssertNoConcurrentBandOverlap(result);
    }

    private static void AssertNoConcurrentBandOverlap(ScheduleResult result)
    {
        for (int i = 0; i < result.Slots.Length; i++)
        {
            for (int j = i + 1; j < result.Slots.Length; j++)
            {
                var a = result.Slots[i];
                var b = result.Slots[j];
                if (!(a.RtStart < b.RtStop && b.RtStart < a.RtStop)) continue;
                double aLo = a.MzMin - Prm / 2.0, aHi = a.MzMax + Prm / 2.0;
                double bLo = b.MzMin - Prm / 2.0, bHi = b.MzMax + Prm / 2.0;
                Assert.False(aLo < bHi && bLo < aHi,
                    $"Concurrent slots overlap: [{aLo:F2}, {aHi:F2}] and [{bLo:F2}, {bHi:F2}]");
            }
        }
    }

    [Fact]
    public void Scheduler_SlotAwarePacking_ConcurrentBandsNeverOverlap()
    {
        // Eight co-eluting precursors 0.5 Th apart, with two clash pairs seeded
        // (p2 clashes p0 on 1000, p4 clashes p1 on 1400) to force window splits.
        var cands = new[]
        {
            Make("p0", 600.0, 5.0, 5.4, "G0", 8e6, new[] { 1000.0, 1100.0, 1200.0, 1300.0 }),
            Make("p1", 600.5, 5.0, 5.4, "G1", 7e6, new[] { 1400.0, 1500.0, 1600.0, 1700.0 }),
            Make("p2", 601.0, 5.0, 5.4, "G2", 6e6, new[] { 1000.0, 1800.0, 1900.0, 2000.0 }),
            Make("p3", 601.5, 5.0, 5.4, "G3", 5e6, new[] { 2100.0, 2200.0, 2300.0, 2400.0 }),
            Make("p4", 602.0, 5.0, 5.4, "G4", 4e6, new[] { 1400.0, 2500.0, 2600.0, 2700.0 }),
            Make("p5", 602.5, 5.0, 5.4, "G5", 3e6, new[] { 2800.0, 2900.0, 3000.0, 3100.0 }),
            Make("p6", 603.0, 5.0, 5.4, "G6", 2e6, new[] { 3200.0, 3300.0, 3400.0, 3500.0 }),
            Make("p7", 603.5, 5.0, 5.4, "G7", 1e6, new[] { 3600.0, 3700.0, 3800.0, 3900.0 }),
        };

        var result = Scheduler.Run(cands, Params());

        for (int i = 0; i < result.Slots.Length; i++)
        {
            for (int j = i + 1; j < result.Slots.Length; j++)
            {
                var a = result.Slots[i];
                var b = result.Slots[j];
                bool concurrent = a.RtStart < b.RtStop && b.RtStart < a.RtStop;
                if (!concurrent) continue;

                double aLo = a.MzMin - Prm / 2.0, aHi = a.MzMax + Prm / 2.0;
                double bLo = b.MzMin - Prm / 2.0, bHi = b.MzMax + Prm / 2.0;
                Assert.False(
                    aLo < bHi && bLo < aHi,
                    $"Concurrent slots overlap in m/z: [{aLo:F2}, {aHi:F2}] and [{bLo:F2}, {bHi:F2}]");
            }
        }
    }
}
