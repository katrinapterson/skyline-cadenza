#nullable enable

using SkylineCadenza.Core.Parsimony;
using SkylineCadenza.Core.Scheduling;

namespace SkylineCadenza.Core.Ingest;

/// <summary>
/// Joins the three ingest products - DIA-NN rows, the parsimony assignment,
/// and the Carafe top-4 fragment table - into the
/// <see cref="Candidate"/> list that <see cref="Scheduler.Run"/> consumes.
/// </summary>
public static class CandidateBuilder
{
    /// <summary>
    /// Build the candidate list. Rows whose peptide isn't assigned to any
    /// parsimonious group are dropped (they don't belong to a coverable
    /// protein). Rows whose Carafe key has no fragment match get an empty
    /// fragment array - the scheduler treats them as solo-slot only.
    /// </summary>
    public static List<Candidate> Build(
        IEnumerable<DiannRow> rows,
        IReadOnlyDictionary<string, (string Group, string Type)> parsimony,
        IReadOnlyDictionary<FragmentKey, FragmentIon[]> fragmentsByKey)
    {
        var candidates = new List<Candidate>();
        foreach (var row in rows)
        {
            if (!parsimony.TryGetValue(row.StrippedSequence, out var assignment))
                continue;

            var carafeKey = new FragmentKey(
                CarafeKey.FromDiann(row.ModifiedSequence),
                row.PrecursorCharge);
            // Confirmed peptides schedule on their MEASURED (DIA-NN Fr.N.Id)
            // fragments, so the empirical schedule is stable regardless of
            // whether a Carafe library is loaded - loading Carafe must not move
            // the confirmed schedule. Carafe's PREDICTED fragments are used only
            // as a fallback for rows that carry no measured fragments (and for
            // predicted peptides, via BuildFromCarafe).
            FragmentIon[] fragments =
                row.Fragments is { Length: > 0 } measured
                    ? measured
                    : fragmentsByKey.TryGetValue(carafeKey, out var ions)
                        ? ions
                        : Array.Empty<FragmentIon>();

            candidates.Add(new Candidate
            {
                PrecursorId = row.PrecursorId,
                StrippedSequence = row.StrippedSequence,
                ModifiedSequence = row.ModifiedSequence,
                PrecursorCharge = row.PrecursorCharge,
                PrecursorMz = row.PrecursorMz,
                RtStart = row.RtStart,
                RtStop = row.RtStop,
                RtApex = row.RtApex,
                PrecursorQuantity = row.PrecursorQuantity,
                Proteotypic = row.Proteotypic,
                QValue = row.QValue,
                ProteinQValue = row.ProteinQValue,
                ProteinGroup = assignment.Group,
                PeptideType = assignment.Type,
                Fragments = fragments,
                Top4Fragments = Candidate.DeriveTopMz(fragments, 4),
                Run = row.Run,
            });
        }
        return candidates;
    }

    /// <summary>
    /// Build candidates from the Carafe library when no DIA-NN report is
    /// available. Peak boundaries are synthesised by padding
    /// <c>TrRecalibrated</c> by <paramref name="peakHalfWidthMin"/> on each
    /// side (the SCHEDULER then layers its 30/15 s firing pad on top).
    /// </summary>
    /// <param name="predicted">
    /// Optional predicted-ranking overlay (Bromo / Carafe detectability). When
    /// supplied, each candidate's <see cref="Candidate.PredictedPriority"/> is
    /// resolved by (peptide, charge) so the scheduler can order peptides within
    /// a protein with <see cref="PeptidePriority.PredictedRanking"/>.
    /// </param>
    /// <param name="detectability">
    /// Optional Carafe detectability overlay. Used as a candidate FILTER
    /// (precursors below <paramref name="minDetectability"/> are dropped;
    /// unknown is kept) AND as the within-protein ranking fallback: when no
    /// <paramref name="predicted"/> (Bromo) overlay is supplied, detectability
    /// populates <see cref="Candidate.PredictedPriority"/> so the scheduler can
    /// order peptides by detectability instead of the PrecursorQuantity stand-in.
    /// </param>
    /// <param name="minDetectability">
    /// Detectability cutoff applied only when <paramref name="detectability"/>
    /// is supplied. Empirically ~0.5 removes peptides that are ~90% unlikely to
    /// be observed while retaining most real detections.
    /// </param>
    public static List<Candidate> BuildFromCarafe(
        IEnumerable<CarafePrecursor> precursors,
        IReadOnlyDictionary<string, (string Group, string Type)> parsimony,
        double peakHalfWidthMin = 0.10,
        PredictedRankingReader.Result? predicted = null,
        PredictedRankingReader.Result? detectability = null,
        double minDetectability = 0.0,
        CandidateSource source = CandidateSource.Experimental)
    {
        bool filterDet = detectability is not null && minDetectability > 0.0;
        var candidates = new List<Candidate>();
        foreach (var p in precursors)
        {
            if (!parsimony.TryGetValue(p.StrippedPeptide, out var assignment))
                continue;
            if (filterDet)
            {
                double? dv = detectability!.Lookup(p.ModifiedPeptide, p.StrippedPeptide, p.PrecursorCharge);
                if (dv.HasValue && dv.Value < minDetectability) continue; // unknown (null) is kept
            }
            string precursorId = $"{p.StrippedPeptide}+{p.PrecursorCharge}";
            candidates.Add(new Candidate
            {
                PrecursorId = precursorId,
                StrippedSequence = p.StrippedPeptide,
                ModifiedSequence = p.ModifiedPeptide,
                PrecursorCharge = p.PrecursorCharge,
                PrecursorMz = p.PrecursorMz,
                RtStart = p.TrRecalibrated - peakHalfWidthMin,
                RtStop = p.TrRecalibrated + peakHalfWidthMin,
                RtApex = p.TrRecalibrated,
                PrecursorQuantity = p.PredictedIntensity,
                // Within-protein priority source, in order: Bromo ranking, then
                // Carafe detectability (reusing the filter overlay), else null
                // (scheduler falls back to the PrecursorQuantity stand-in).
                PredictedPriority =
                    predicted?.Lookup(p.ModifiedPeptide, p.StrippedPeptide, p.PrecursorCharge)
                    ?? detectability?.Lookup(p.ModifiedPeptide, p.StrippedPeptide, p.PrecursorCharge),
                QValue = 0.0,
                ProteinQValue = double.NaN,
                Proteotypic = p.ProteinAccessions.Length == 1 ? 1 : 0,
                ProteinGroup = assignment.Group,
                PeptideType = assignment.Type,
                Fragments = p.Fragments,
                Top4Fragments = Candidate.DeriveTopMz(p.Fragments, 4),
                Run = "Carafe (predicted)",
                Source = source,
            });
        }
        return candidates;
    }

    /// <summary>
    /// Combined DIA-NN + Carafe "backfill" pool. Experimental candidates come
    /// from the DIA-NN report (tagged <see cref="CandidateSource.Experimental"/>
    /// via <see cref="Build"/>); predicted candidates come from the Carafe
    /// library (tagged <see cref="CandidateSource.Predicted"/>). Predicted
    /// precursors that duplicate a confirmed one on (stripped sequence, charge)
    /// are dropped, so the scheduler never sees the same target twice. The
    /// scheduler then schedules the experimental tier to exhaustion before the
    /// predicted tier competes for leftover budget.
    /// </summary>
    /// <param name="parsimony">
    /// Peptide -> (group, type) assignment covering both tiers. To keep the
    /// confirmed schedule locked, build the DIA-NN peptides' entries from a
    /// DIA-NN-only parsimony run and route the Carafe peptides onto those groups
    /// with <see cref="MapPredictedToExperimentalGroups"/> (rather than a single
    /// parsimony run over the union, which would re-group the confirmed
    /// evidence). Confirmed entries must win for any peptide seen by both tiers.
    /// </param>
    public static List<Candidate> BuildCombined(
        IEnumerable<DiannRow> rows,
        IReadOnlyDictionary<FragmentKey, FragmentIon[]> fragmentsByKey,
        IEnumerable<CarafePrecursor> carafePrecursors,
        IReadOnlyDictionary<string, (string Group, string Type)> parsimony,
        double peakHalfWidthMin = 0.10,
        PredictedRankingReader.Result? predicted = null,
        PredictedRankingReader.Result? detectability = null,
        double minDetectability = 0.0)
    {
        var experimental = Build(rows, parsimony, fragmentsByKey);
        var confirmed = new HashSet<(string, int)>(
            experimental.Select(c => (c.StrippedSequence, c.PrecursorCharge)));

        var predictedCandidates = BuildFromCarafe(
            carafePrecursors, parsimony, peakHalfWidthMin,
            predicted, detectability, minDetectability,
            source: CandidateSource.Predicted);

        var combined = new List<Candidate>(experimental.Count + predictedCandidates.Count);
        combined.AddRange(experimental);
        foreach (var c in predictedCandidates)
            if (!confirmed.Contains((c.StrippedSequence, c.PrecursorCharge)))
                combined.Add(c);
        return combined;
    }

    /// <summary>
    /// Build the peptide -> (group, type) map for the predicted (Carafe) tier
    /// WITHOUT re-running parsimony over the combined evidence, so the confirmed
    /// (DIA-NN) grouping is never perturbed.
    ///
    /// Each Carafe peptide is routed by its protein accessions:
    /// <list type="bullet">
    /// <item>If any accession was seen by DIA-NN, the peptide joins that
    /// experimental canonical group (lexicographically smallest when its
    /// accessions span several confirmed groups), so it counts toward the same
    /// per-protein cap as the confirmed peptides there. It is typed
    /// <c>unique</c> when it maps to a single confirmed group, else
    /// <c>razor</c>.</item>
    /// <item>Peptides whose proteins are ALL unseen by DIA-NN (target-list-only
    /// proteins) are grouped among themselves by a parsimony run over just those
    /// peptides, so predicted-only proteins still get sensible groups.</item>
    /// </list>
    /// </summary>
    public static Dictionary<string, (string Group, string Type)> MapPredictedToExperimentalGroups(
        IReadOnlyDictionary<string, IReadOnlyList<string>> carafePeptideProteins,
        IReadOnlyDictionary<string, string> experimentalAccessionToGroup)
    {
        var result = new Dictionary<string, (string Group, string Type)>();
        var unseen = new Dictionary<string, IReadOnlyList<string>>();
        foreach (var (peptide, accessions) in carafePeptideProteins)
        {
            var seenGroups = new SortedSet<string>(StringComparer.Ordinal);
            foreach (var acc in accessions)
                if (experimentalAccessionToGroup.TryGetValue(acc, out var g))
                    seenGroups.Add(g);

            if (seenGroups.Count > 0)
                result[peptide] = (seenGroups.Min!, seenGroups.Count == 1 ? "unique" : "razor");
            else
                unseen[peptide] = accessions;
        }

        if (unseen.Count > 0)
            foreach (var (peptide, groupType) in Parsimony.ParsimonyEngine.Assign(unseen))
                result[peptide] = groupType;

        return result;
    }

    /// <summary>
    /// Union two peptide->proteins maps for a single parsimony run over the
    /// combined evidence. Protein sets for a shared peptide are merged.
    /// </summary>
    public static Dictionary<string, IReadOnlyList<string>> MergePeptideProteinMaps(
        IReadOnlyDictionary<string, IReadOnlyList<string>> a,
        IReadOnlyDictionary<string, IReadOnlyList<string>> b)
    {
        var merged = new Dictionary<string, HashSet<string>>();
        void Add(IReadOnlyDictionary<string, IReadOnlyList<string>> src)
        {
            foreach (var (peptide, prots) in src)
            {
                if (!merged.TryGetValue(peptide, out var set))
                {
                    set = new HashSet<string>();
                    merged[peptide] = set;
                }
                foreach (var p in prots) set.Add(p);
            }
        }
        Add(a);
        Add(b);
        return merged.ToDictionary(
            kv => kv.Key,
            kv => (IReadOnlyList<string>)kv.Value.ToList());
    }

    /// <summary>
    /// Build a peptide -> proteins map from Carafe precursors for parsimony.
    /// Multiple Carafe rows for the same stripped peptide get their
    /// accession sets unioned.
    /// </summary>
    public static Dictionary<string, IReadOnlyList<string>> CarafePeptideProteinMap(
        IEnumerable<CarafePrecursor> precursors)
    {
        var byPeptide = new Dictionary<string, HashSet<string>>();
        foreach (var p in precursors)
        {
            if (!byPeptide.TryGetValue(p.StrippedPeptide, out var set))
            {
                set = new HashSet<string>();
                byPeptide[p.StrippedPeptide] = set;
            }
            foreach (var acc in p.ProteinAccessions) set.Add(acc);
        }
        return byPeptide.ToDictionary(
            kv => kv.Key,
            kv => (IReadOnlyList<string>)kv.Value.ToList());
    }
}
