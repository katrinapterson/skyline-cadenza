using SkylineCadenza.Core.Scheduling;

namespace SkylineCadenza.Core.Ingest;

/// <summary>
/// What a DIA-NN report offers the scheduler, per protein group, in the
/// names the scheduler uses. Lets another program (AutoTargeting, and Claude
/// through it) choose target proteins and forced peptides that will match:
/// the groups here are the ones <see cref="SchedulingParameters.TargetProteins"/>
/// is compared against, after the same q-value filter and parsimony.
/// </summary>
public static class ReportCatalog
{
    /// <summary>One precursor the scheduler could place.</summary>
    public sealed record PrecursorEntry(
        string ModifiedSequence,
        string StrippedSequence,
        int Charge,
        double Mz,
        double RtApex,
        double RtStart,
        double RtStop,
        double Quantity,
        double QValue,
        string PeptideType);

    /// <summary>One protein group and its precursors, strongest first.</summary>
    public sealed record ProteinEntry(
        string Group,
        IReadOnlyList<string> Genes,
        int Peptides,
        int UniquePeptides,
        double SummedQuantity,
        double? ProteinQValue,
        IReadOnlyList<PrecursorEntry> Precursors);

    /// <summary>
    /// Every protein group in <paramref name="candidates"/>, by summed
    /// precursor quantity (the app's default protein ranking), largest first.
    /// </summary>
    public static List<ProteinEntry> Build(
        IReadOnlyList<Candidate> candidates,
        IReadOnlyDictionary<string, IReadOnlyList<string>> groupGenes)
    {
        return candidates
            .GroupBy(c => c.ProteinGroup, StringComparer.Ordinal)
            .Select(g => new ProteinEntry(
                g.Key,
                groupGenes.TryGetValue(g.Key, out var genes) ? genes : Array.Empty<string>(),
                g.Select(c => c.StrippedSequence).Distinct(StringComparer.Ordinal).Count(),
                g.Where(c => c.PeptideType == "unique").Select(c => c.StrippedSequence).Distinct(StringComparer.Ordinal).Count(),
                g.Sum(c => c.PrecursorQuantity),
                // Null, not NaN, when the report has no PG.Q.Value: JSON has no NaN.
                g.Select(c => c.ProteinQValue).Where(q => !double.IsNaN(q)).Select(q => (double?)q).Min(),
                g.OrderByDescending(c => c.PrecursorQuantity)
                    .Select(c => new PrecursorEntry(
                        c.ModifiedSequence, c.StrippedSequence, c.PrecursorCharge, c.PrecursorMz,
                        c.RtApex, c.RtStart, c.RtStop, c.PrecursorQuantity, c.QValue, c.PeptideType))
                    .ToList()))
            .OrderByDescending(p => p.SummedQuantity)
            .ThenBy(p => p.Group, StringComparer.Ordinal)
            .ToList();
    }

    /// <summary>
    /// Gene names for each protein group in <paramref name="candidates"/>.
    /// DIA-NN lists <c>Protein.Ids</c> and <c>Genes</c> semicolon-separated;
    /// when the two lists are the same length they are paired by position,
    /// and when a row names a single gene it applies to all of its proteins.
    /// Other rows are ambiguous and skipped. A group is one canonical
    /// accession, so its genes are that accession's.
    /// </summary>
    public static Dictionary<string, IReadOnlyList<string>> GroupGenes(
        IEnumerable<DiannRow> rows, IEnumerable<Candidate> candidates)
    {
        var byAccession = new Dictionary<string, SortedSet<string>>(StringComparer.Ordinal);
        foreach (var row in rows)
        {
            var accs = row.ProteinIds.Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
            var genes = row.Genes.Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
            if (genes.Length == 0) continue;
            for (int i = 0; i < accs.Length; i++)
            {
                var gene = genes.Length == accs.Length ? genes[i] : genes.Length == 1 ? genes[0] : null;
                if (gene is null) continue;
                if (!byAccession.TryGetValue(accs[i], out var set))
                    byAccession[accs[i]] = set = new SortedSet<string>(StringComparer.Ordinal);
                set.Add(gene);
            }
        }

        var result = new Dictionary<string, IReadOnlyList<string>>(StringComparer.Ordinal);
        foreach (var group in candidates.Select(c => c.ProteinGroup).Distinct(StringComparer.Ordinal))
        {
            if (byAccession.TryGetValue(group, out var set)
                || byAccession.TryGetValue(ProteinListParser.ExtractAccession(group), out set))
                result[group] = set.ToList();
        }
        return result;
    }
}
