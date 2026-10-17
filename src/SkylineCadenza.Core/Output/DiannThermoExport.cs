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

        /// <summary>
        /// Each target-list entry, in order, with the accessions it resolved
        /// to (a gene can name several). Empty when no list was given.
        /// </summary>
        public IReadOnlyList<(string Entry, IReadOnlyList<string> Accessions)> TargetEntries { get; init; } =
            Array.Empty<(string, IReadOnlyList<string>)>();

        /// <summary>Gene names per protein group, for summaries.</summary>
        public IReadOnlyDictionary<string, IReadOnlyList<string>> GroupGenes { get; init; } =
            new Dictionary<string, IReadOnlyList<string>>();
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
        var (rows, parsimony, candidates) = await LoadPoolAsync(reportPath, parameters.QValueCutoff, cancellationToken);
        var geneToAccession = ProteinListParser.BuildGeneToAccession(rows);
        var knownAccessions = parsimony.Values.Select(v => v.Group).ToHashSet();

        var p = parameters;
        var entries = new List<(string, IReadOnlyList<string>)>();
        if (!string.IsNullOrWhiteSpace(targetProteinsText))
        {
            var tokens = ProteinListParser.ParseText(targetProteinsText);
            var resolved = ProteinListParser.ResolveGenes(tokens, geneToAccession, knownAccessions);
            p = p with { TargetProteins = resolved.ToHashSet() };
            // Per entry, so a caller can say which ones matched nothing.
            foreach (var token in tokens)
                entries.Add((token, ProteinListParser.ResolveGenes(new[] { token }, geneToAccession, knownAccessions)));
        }

        var schedule = Scheduler.Run(candidates, p, cancellationToken);
        var csv = ThermoCsvWriter.Build(candidates, schedule, p);
        return new Result(candidates, schedule, p, csv)
        {
            TargetEntries = entries,
            GroupGenes = ReportCatalog.GroupGenes(rows, candidates),
        };
    }

    /// <summary>
    /// The protein groups and precursors the scheduler would choose from in
    /// <paramref name="reportPath"/>: the same ingest as <see cref="RunAsync"/>,
    /// without scheduling.
    /// </summary>
    public static async Task<List<ReportCatalog.ProteinEntry>> DescribeAsync(
        string reportPath, double qValueCutoff, CancellationToken cancellationToken = default)
    {
        var (rows, _, candidates) = await LoadPoolAsync(reportPath, qValueCutoff, cancellationToken);
        return ReportCatalog.Build(candidates, ReportCatalog.GroupGenes(rows, candidates));
    }

    private static async Task<(List<DiannRow> Rows, Dictionary<string, (string Group, string Type)> Parsimony, List<Candidate> Candidates)>
        LoadPoolAsync(string reportPath, double qValueCutoff, CancellationToken cancellationToken)
    {
        var rows = await DiannParquetReader.LoadAsync(reportPath, qValueCutoff, cancellationToken);
        if (rows.Count == 0)
            throw new InvalidDataException(
                $"No precursors in {reportPath} passed the q-value cutoff {qValueCutoff}.");

        var parsimony = ParsimonyEngine.Assign(DiannParquetReader.BuildPeptideProteinMap(rows));
        var candidates = CandidateBuilder.Build(rows, parsimony, new Dictionary<FragmentKey, FragmentIon[]>());
        return (rows, parsimony, candidates);
    }
}
