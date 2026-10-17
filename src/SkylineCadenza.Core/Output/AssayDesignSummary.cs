using SkylineCadenza.Core.Ingest;
using SkylineCadenza.Core.Scheduling;

namespace SkylineCadenza.Core.Output;

/// <summary>
/// What an <see cref="DiannThermoExport"/> run scheduled, for a caller that
/// cannot look at the app's plots: how much of the report made it in, how
/// full the cycle is, which target-list entries and forced peptides were
/// placed, and the peptides chosen per protein. Written next to the CSV by
/// <c>SkylineCadenzaCli export-thermo ... SUMMARY.json</c>.
/// </summary>
public sealed record AssayDesignSummary
{
    public required int PoolPrecursors { get; init; }
    public required int PoolProteinGroups { get; init; }
    public required int ScheduledPrecursors { get; init; }
    public required int ScheduledPeptides { get; init; }

    /// <summary>Precursors that ride inside a scheduled window without being targeted themselves.</summary>
    public required int Passengers { get; init; }

    public required int ProteinGroupsCovered { get; init; }
    public required int Windows { get; init; }
    public required int CsvRows { get; init; }
    public required string Mode { get; init; }
    public required int CycleBudget { get; init; }

    /// <summary>Most windows open at one retention time; compare with <see cref="CycleBudget"/>.</summary>
    public required int PeakConcurrentWindows { get; init; }

    public required IReadOnlyList<TargetOutcome> Targets { get; init; }
    public required IReadOnlyList<ForcedOutcome> Forced { get; init; }
    public required IReadOnlyList<ProteinOutcome> Proteins { get; init; }

    /// <param name="Entry">The target-list entry as given (accession or gene).</param>
    /// <param name="ProteinGroups">The report's protein groups it matched; empty when it matched none.</param>
    /// <param name="ScheduledPrecursors">Precursors of those groups in the schedule.</param>
    public sealed record TargetOutcome(string Entry, IReadOnlyList<string> ProteinGroups, int ScheduledPrecursors);

    /// <param name="Status">Scheduled, NotInPool (no identification) or Rejected (failed the scheduling rules).</param>
    public sealed record ForcedOutcome(string Entry, string Status, string Precursor, bool UsedFallbackCharge);

    public sealed record ProteinOutcome(
        string Group, IReadOnlyList<string> Genes, bool IsTarget, IReadOnlyList<string> Precursors);

    public static AssayDesignSummary Build(DiannThermoExport.Result export)
    {
        var candidates = export.Candidates;
        var schedule = export.Schedule;
        var p = export.Parameters;
        var scheduled = schedule.ScheduledIndices.Select(i => candidates[i]).ToList();

        bool Matches(IEnumerable<string> accessions, string group)
        {
            string accession = ProteinListParser.ExtractAccession(group);
            return accessions.Any(a => a == group || a == accession);
        }

        var poolGroups = candidates.Select(c => c.ProteinGroup).Distinct(StringComparer.Ordinal).ToList();
        var targets = export.TargetEntries
            .Select(t =>
            {
                var groups = poolGroups.Where(g => Matches(t.Accessions, g)).ToList();
                return new TargetOutcome(t.Entry, groups, scheduled.Count(c => groups.Contains(c.ProteinGroup)));
            })
            .ToList();

        var proteins = scheduled
            .GroupBy(c => c.ProteinGroup, StringComparer.Ordinal)
            .Select(g => new ProteinOutcome(
                g.Key,
                export.GroupGenes.TryGetValue(g.Key, out var genes) ? genes : Array.Empty<string>(),
                p.TargetProteins.Count > 0 && Matches(p.TargetProteins, g.Key),
                g.OrderByDescending(c => c.PrecursorQuantity).Select(c => $"{c.ModifiedSequence}+{c.PrecursorCharge}").ToList()))
            .OrderByDescending(x => x.IsTarget)
            .ThenBy(x => x.Group, StringComparer.Ordinal)
            .ToList();

        return new AssayDesignSummary
        {
            PoolPrecursors = candidates.Count,
            PoolProteinGroups = poolGroups.Count,
            ScheduledPrecursors = scheduled.Count,
            ScheduledPeptides = scheduled.Select(c => c.StrippedSequence).Distinct(StringComparer.Ordinal).Count(),
            Passengers = schedule.PassengerIndices.Length,
            ProteinGroupsCovered = schedule.ProteinGroupsCovered,
            Windows = schedule.Slots.Length,
            CsvRows = export.RowCount,
            Mode = p.Mode.ToString(),
            CycleBudget = p.CycleBudget,
            PeakConcurrentWindows = schedule.SlotCountCurve.DefaultIfEmpty(0).Max(),
            Targets = targets,
            Forced = schedule.ForcedOutcomes
                .Select(o => new ForcedOutcome(
                    o.Spec.Text,
                    o.Status.ToString(),
                    o.CandidateIndex >= 0 ? $"{candidates[o.CandidateIndex].ModifiedSequence}+{candidates[o.CandidateIndex].PrecursorCharge}" : null,
                    o.UsedFallbackCharge))
                .ToList(),
            Proteins = proteins,
        };
    }
}
