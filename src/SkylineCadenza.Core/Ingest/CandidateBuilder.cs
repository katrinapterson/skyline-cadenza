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
            // Prefer the user-supplied Carafe fragments if present.
            // Otherwise fall back to DIA-NN's per-row Fr.N.Id fragments
            // (real measured m/z), which gives a useful interference check
            // without requiring the user to load Carafe.
            FragmentIon[] fragments =
                fragmentsByKey.TryGetValue(carafeKey, out var ions)
                    ? ions
                    : row.Fragments ?? Array.Empty<FragmentIon>();

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
    /// Parsimony assignment built over the UNION of the DIA-NN and Carafe
    /// peptide->protein maps (see <see cref="MergePeptideProteinMaps"/>) so a
    /// predicted peptide lands in the same protein group as the confirmed
    /// evidence and counts toward that group's per-protein cap.
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
