using SkylineCadenza.Core.Ingest;
using SkylineCadenza.Core.Parsimony;
using SkylineCadenza.Core.Scheduling;
using SkylineCadenza.Core.SkylineRpc;

namespace SkylineCadenza.Core.Output;

/// <summary>
/// Headless DIA-NN report -> Thermo scheduled inclusion CSV. Runs the same
/// Core calls, in the same order, as the app's Empirical (DIA-NN report) Load
/// followed by "Export Thermo CSV...": <see cref="DiannParquetReader.LoadAsync"/>
/// at the q-value cutoff, <see cref="ParsimonyEngine.Assign"/>,
/// <see cref="CandidateBuilder.Build"/> on the report's own fragments,
/// target-list gene resolution, <see cref="Scheduler.Run"/>, and
/// <see cref="ThermoCsvWriter.Build"/>. Used by <c>SkylineCadenzaCli
/// export-thermo</c> so another app (AutoTargeting) gets the CSV the GUI
/// would write with the same settings.
/// </summary>
public static class DiannThermoExport
{
    /// <summary>
    /// The app's defaults for a fresh window. These differ from a bare
    /// <see cref="SchedulingParameters"/> in two places: the app starts with
    /// <see cref="SchedulingParameters.MinTargetIonShare"/> at 0.015 and always
    /// runs with <see cref="SchedulingParameters.BackgroundFragmentClash"/> on.
    /// Every other knob uses the record's own default.
    /// </summary>
    public static SchedulingParameters AppDefaults { get; } = new()
    {
        MinTargetIonShare = 0.015,
        BackgroundFragmentClash = true,
    };

    /// <summary>Everything the export produced, for callers that want more than the CSV.</summary>
    public sealed record Result(
        IReadOnlyList<Candidate> Candidates,
        ScheduleResult Schedule,
        SchedulingParameters Parameters,
        string Csv)
    {
        /// <summary>Data rows in <see cref="Csv"/> (excludes the header).</summary>
        public int RowCount => Csv.Split('\n', StringSplitOptions.RemoveEmptyEntries).Length - 1;
    }

    /// <summary>
    /// Ingest <paramref name="reportPath"/>, schedule with
    /// <paramref name="parameters"/>, and build the Thermo CSV.
    /// </summary>
    /// <param name="targetProteinsText">
    /// Optional target protein list in the app's box syntax (accessions or
    /// gene names, FASTA headers, one per line or separated by whitespace /
    /// commas). Resolved against the report exactly as the app does and
    /// written into <see cref="SchedulingParameters.TargetProteins"/>.
    /// </param>
    public static async Task<Result> RunAsync(
        string reportPath,
        SchedulingParameters parameters,
        string? targetProteinsText = null,
        CancellationToken cancellationToken = default)
    {
        var rows = await DiannParquetReader.LoadAsync(reportPath, parameters.QValueCutoff, cancellationToken);
        if (rows.Count == 0)
            throw new InvalidDataException(
                $"No precursors in {reportPath} passed the q-value cutoff {parameters.QValueCutoff}.");

        var geneToAccession = ProteinListParser.BuildGeneToAccession(rows);
        var parsimony = ParsimonyEngine.Assign(DiannParquetReader.BuildPeptideProteinMap(rows));
        var knownAccessions = parsimony.Values.Select(v => v.Group).ToHashSet();
        var candidates = CandidateBuilder.Build(rows, parsimony, new Dictionary<FragmentKey, FragmentIon[]>());

        var p = parameters;
        if (!string.IsNullOrWhiteSpace(targetProteinsText))
        {
            var tokens = ProteinListParser.ParseText(targetProteinsText);
            var resolved = ProteinListParser.ResolveGenes(tokens, geneToAccession, knownAccessions);
            p = p with { TargetProteins = resolved.ToHashSet() };
        }

        var schedule = Scheduler.Run(candidates, p, cancellationToken);
        var csv = ThermoCsvWriter.Build(candidates, schedule, p);
        return new Result(candidates, schedule, p, csv);
    }
}
