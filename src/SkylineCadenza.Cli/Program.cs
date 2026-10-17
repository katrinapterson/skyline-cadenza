using System.Text.Json;
using System.Text.Json.Serialization;
using SkylineCadenza.Core.Ingest;
using SkylineCadenza.Core.Output;
using SkylineCadenza.Core.Parsimony;
using SkylineCadenza.Core.Scheduling;

namespace SkylineCadenza.Cli;

/// <summary>
/// Headless companion to the WPF app. Reads a candidate pool + a
/// <see cref="SchedulingParameters"/> from JSON, runs
/// <see cref="Scheduler.Run"/>, and writes the schedule back out as
/// JSON. Used by the algorithm-comparison Jupyter notebook so that the
/// realistic plots reflect Cadenza's real C# scheduler bit-for-bit
/// rather than a Python re-implementation. <c>export-thermo</c> writes the
/// app's Thermo scheduled inclusion CSV from a DIA-NN report, so another app
/// (AutoTargeting) can design an assay without opening the GUI.
/// </summary>
public static class Program
{
    private static readonly JsonSerializerOptions JsonOpts = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
        WriteIndented = true,
        Converters = { new JsonStringEnumConverter() },
    };

    /// <summary>
    /// Entry point. Usage:
    ///   SkylineCadenzaCli schedule INPUT.json OUTPUT.json
    ///   SkylineCadenzaCli schedule - -        (stdin / stdout)
    ///   SkylineCadenzaCli schedule-diann REPORT.parquet PARAMS.json OUTPUT.json
    ///   SkylineCadenzaCli export-thermo REPORT.parquet PARAMS.json|- OUTPUT.csv [SUMMARY.json] [--skyline OUT.sky] [--skylinecmd PATH]
    ///   SkylineCadenzaCli describe-diann REPORT.parquet PARAMS.json|- OUTPUT.json
    /// </summary>
    public static async Task<int> Main(string[] args)
    {
        try
        {
            if (args.Length < 1)
            {
                PrintUsage();
                return 64; // EX_USAGE
            }

            return args[0] switch
            {
                "schedule" => await RunScheduleAsync(args.Skip(1).ToArray()),
                "schedule-diann" => await RunScheduleDiannAsync(args.Skip(1).ToArray()),
                "export-thermo" => await RunExportThermoAsync(args.Skip(1).ToArray()),
                "describe-diann" => await RunDescribeDiannAsync(args.Skip(1).ToArray()),
                "--help" or "-h" or "help" => Help(),
                _ => UnknownVerb(args[0]),
            };
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine($"error: {ex.GetType().Name}: {ex.Message}");
            return 1;
        }
    }

    private static async Task<int> RunScheduleAsync(string[] rest)
    {
        if (rest.Length != 2)
        {
            Console.Error.WriteLine("usage: SkylineCadenzaCli schedule INPUT OUTPUT  (use '-' for stdin/stdout)");
            return 64;
        }

        ScheduleRequest request = rest[0] == "-"
            ? await JsonSerializer.DeserializeAsync<ScheduleRequest>(Console.OpenStandardInput(), JsonOpts)
              ?? throw new InvalidDataException("Empty stdin")
            : JsonSerializer.Deserialize<ScheduleRequest>(File.ReadAllText(rest[0]), JsonOpts)
              ?? throw new InvalidDataException($"Empty file: {rest[0]}");

        var candidates = request.Candidates.Select(JsonToCandidate).ToList();
        var parameters = ToSchedulingParameters(request.Parameters, new SchedulingParameters());
        var result = Scheduler.Run(candidates, parameters);
        var response = ToResponse(result);

        if (rest[1] == "-")
        {
            await JsonSerializer.SerializeAsync(Console.OpenStandardOutput(), response, JsonOpts);
            Console.Out.WriteLine();
        }
        else
        {
            await File.WriteAllTextAsync(rest[1], JsonSerializer.Serialize(response, JsonOpts));
        }
        return 0;
    }

    /// <summary>
    /// Ingest a DIA-NN report.parquet through the same path the app's
    /// DIA-NN-only Load uses (reader at the parameters' q-value cutoff,
    /// parsimony, <see cref="CandidateBuilder.Build"/> on the report's own
    /// fragments), schedule it, and write one row per scheduled precursor.
    /// Lets offline analyses compare parameter settings on a real report
    /// without re-implementing ingest.
    /// </summary>
    private static async Task<int> RunScheduleDiannAsync(string[] rest)
    {
        if (rest.Length != 3)
        {
            Console.Error.WriteLine("usage: SkylineCadenzaCli schedule-diann REPORT.parquet PARAMS.json OUTPUT.json");
            return 64;
        }

        var paramDto = JsonSerializer.Deserialize<SchedulingParametersDto>(File.ReadAllText(rest[1]), JsonOpts);
        var parameters = ToSchedulingParameters(paramDto, new SchedulingParameters());
        var rows = await DiannParquetReader.LoadAsync(rest[0], parameters.QValueCutoff);
        var parsimony = ParsimonyEngine.Assign(DiannParquetReader.BuildPeptideProteinMap(rows));
        var candidates = CandidateBuilder.Build(rows, parsimony, new Dictionary<FragmentKey, FragmentIon[]>());
        var result = Scheduler.Run(candidates, parameters);

        var targets = new List<ScheduledTargetDto>(result.ScheduledIndices.Length);
        for (int k = 0; k < result.ScheduledIndices.Length; k++)
        {
            var c = candidates[result.ScheduledIndices[k]];
            var slot = result.Slots[result.ScheduledSlotIds[k]];
            targets.Add(new ScheduledTargetDto
            {
                PrecursorId = c.PrecursorId,
                StrippedSequence = c.StrippedSequence,
                PrecursorCharge = c.PrecursorCharge,
                ProteinGroup = c.ProteinGroup,
                PrecursorMz = c.PrecursorMz,
                PrecursorQuantity = c.PrecursorQuantity,
                SlotId = slot.Id,
                SlotMzMin = slot.MzMin,
                SlotMzMax = slot.MzMax,
                SlotRtStart = slot.RtStart,
                SlotRtStop = slot.RtStop,
                SlotSize = slot.MemberIndices.Count,
            });
        }
        var response = new DiannScheduleResponse
        {
            PoolPrecursors = candidates.Count,
            ProteinGroupsCovered = result.ProteinGroupsCovered,
            Slots = result.Slots.Length,
            Targets = targets.ToArray(),
            Forced = result.ForcedOutcomes.Select(o => new ForcedOutcomeDto
            {
                Entry = o.Spec.Text,
                Status = o.Status.ToString(),
                PrecursorId = o.CandidateIndex >= 0 ? candidates[o.CandidateIndex].PrecursorId : null,
                UsedFallbackCharge = o.UsedFallbackCharge,
            }).ToArray(),
        };
        await File.WriteAllTextAsync(rest[2], JsonSerializer.Serialize(response, JsonOpts));
        Console.Error.WriteLine(
            $"pool {candidates.Count:n0}  scheduled {targets.Count:n0}  proteins {result.ProteinGroupsCovered:n0}  slots {result.Slots.Length:n0}");
        return 0;
    }

    /// <summary>
    /// Ingest a DIA-NN report.parquet and write the Thermo Method Editor
    /// scheduled inclusion list the app's "Export Thermo CSV..." writes, via
    /// <see cref="DiannThermoExport.RunAsync"/> (the app's Empirical Load +
    /// export path). Unset parameters take the APP's defaults
    /// (<see cref="DiannThermoExport.AppDefaults"/>), not the bare Core ones,
    /// so a headless run matches the GUI with untouched settings. PARAMS may
    /// be '-' for all defaults. An optional SUMMARY.json gets an
    /// <see cref="AssayDesignSummary"/>, so a caller can see what was
    /// scheduled without parsing the CSV; it is written even when nothing was
    /// scheduled, to say why.
    /// </summary>
    internal static async Task<int> RunExportThermoAsync(string[] rest)
    {
        if (!TryTakeOptions(rest, out var positional, out var options, "--skyline", "--skylinecmd")
            || positional.Count is not (3 or 4))
        {
            Console.Error.WriteLine("usage: SkylineCadenzaCli export-thermo REPORT.parquet PARAMS.json|- OUTPUT.csv [SUMMARY.json] "
                + "[--skyline OUTPUT.sky] [--skylinecmd SkylineCmd.exe]");
            return 64;
        }
        string reportPath = positional[0], paramsPath = positional[1], outputPath = positional[2];
        string? summaryPath = positional.Count == 4 ? positional[3] : null;
        options.TryGetValue("--skyline", out var skylinePath);
        options.TryGetValue("--skylinecmd", out var skylineCmd);
        if (!File.Exists(reportPath))
        {
            Console.Error.WriteLine($"error: DIA-NN report not found: {reportPath}");
            return 66; // EX_NOINPUT
        }
        if (paramsPath != "-" && !File.Exists(paramsPath))
        {
            Console.Error.WriteLine($"error: parameters file not found: {paramsPath}");
            return 66;
        }

        var (parameters, targetText) = ReadExportParameters(paramsPath == "-" ? null : File.ReadAllText(paramsPath));
        var export = await DiannThermoExport.RunAsync(reportPath, parameters, targetText);
        if (summaryPath is not null)
            await WriteJsonAsync(summaryPath, AssayDesignSummary.Build(export));
        if (export.Schedule.ScheduledIndices.Length == 0)
        {
            Console.Error.WriteLine(
                $"error: nothing scheduled from {export.Candidates.Count:n0} precursors; no CSV written. Check the parameters.");
            return 2;
        }

        string? dir = Path.GetDirectoryName(Path.GetFullPath(outputPath));
        if (!string.IsNullOrEmpty(dir)) Directory.CreateDirectory(dir);
        await File.WriteAllTextAsync(outputPath, export.Csv);

        var r = export.Schedule;
        Console.Error.WriteLine(
            $"pool {export.Candidates.Count:n0}  scheduled {r.ScheduledIndices.Length:n0}  proteins {r.ProteinGroupsCovered:n0}  "
            + $"slots {r.Slots.Length:n0}  rows {export.RowCount:n0}  mode {export.Parameters.Mode}  -> {outputPath}");

        if (skylinePath is not null)
        {
            try
            {
                var document = await SkylineDocumentExport.WriteAsync(export, skylinePath, skylineCmd);
                Console.Error.WriteLine($"skyline {document.Precursors:n0} precursors  -> {document.Document}");
            }
            catch (InvalidOperationException ex)
            {
                // The CSV stands; only the document is missing.
                Console.Error.WriteLine($"error: {ex.Message}");
                return 3;
            }
        }
        return 0;
    }

    /// <summary>
    /// Split <paramref name="args"/> into positional arguments and the named
    /// <paramref name="names"/> (<c>--name value</c> or <c>--name=value</c>).
    /// False on an unknown option or one without a value.
    /// </summary>
    internal static bool TryTakeOptions(
        IReadOnlyList<string> args, out List<string> positional, out Dictionary<string, string> options, params string[] names)
    {
        positional = new List<string>();
        options = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        for (int i = 0; i < args.Count; i++)
        {
            var arg = args[i];
            if (!arg.StartsWith("--", StringComparison.Ordinal))
            {
                positional.Add(arg);
                continue;
            }
            var eq = arg.IndexOf('=');
            var name = eq > 0 ? arg[..eq] : arg;
            if (!names.Contains(name, StringComparer.OrdinalIgnoreCase)) return false;
            if (eq > 0) options[name] = arg[(eq + 1)..];
            else if (i + 1 < args.Count) options[name] = args[++i];
            else return false;
        }
        return true;
    }

    /// <summary>
    /// Write the protein groups and precursors a DIA-NN report offers the
    /// scheduler (<see cref="ReportCatalog"/>), after export-thermo's ingest
    /// at PARAMS' q-value cutoff, so a caller can pick target proteins and
    /// forced peptides by the names export-thermo will match.
    /// </summary>
    internal static async Task<int> RunDescribeDiannAsync(string[] rest)
    {
        if (rest.Length != 3)
        {
            Console.Error.WriteLine("usage: SkylineCadenzaCli describe-diann REPORT.parquet PARAMS.json|- OUTPUT.json");
            return 64;
        }
        string reportPath = rest[0], paramsPath = rest[1], outputPath = rest[2];
        if (!File.Exists(reportPath))
        {
            Console.Error.WriteLine($"error: DIA-NN report not found: {reportPath}");
            return 66;
        }
        if (paramsPath != "-" && !File.Exists(paramsPath))
        {
            Console.Error.WriteLine($"error: parameters file not found: {paramsPath}");
            return 66;
        }

        var (parameters, _) = ReadExportParameters(paramsPath == "-" ? null : File.ReadAllText(paramsPath));
        var proteins = await DiannThermoExport.DescribeAsync(reportPath, parameters.QValueCutoff);
        await WriteJsonAsync(outputPath, new ReportDescription
        {
            QValueCutoff = parameters.QValueCutoff,
            ProteinGroups = proteins.Count,
            Precursors = proteins.Sum(p => p.Precursors.Count),
            Proteins = proteins,
        });
        Console.Error.WriteLine($"proteins {proteins.Count:n0}  precursors {proteins.Sum(p => p.Precursors.Count):n0}  -> {outputPath}");
        return 0;
    }

    private static async Task WriteJsonAsync<T>(string path, T value)
    {
        string? dir = Path.GetDirectoryName(Path.GetFullPath(path));
        if (!string.IsNullOrEmpty(dir)) Directory.CreateDirectory(dir);
        await File.WriteAllTextAsync(path, JsonSerializer.Serialize(value, JsonOpts));
    }

    /// <summary>
    /// Parse export-thermo's PARAMS JSON (null or blank = all defaults) over
    /// the app's defaults. Returns the parameters plus the raw target-protein
    /// text, which is resolved against the report during the export.
    /// </summary>
    internal static (SchedulingParameters Parameters, string? TargetProteinsText) ReadExportParameters(string? json)
    {
        var dto = string.IsNullOrWhiteSpace(json)
            ? null
            : JsonSerializer.Deserialize<SchedulingParametersDto>(json, JsonOpts);
        var p = ToSchedulingParameters(dto, DiannThermoExport.AppDefaults);
        string? targets = dto?.TargetProteins is { Length: > 0 } t ? string.Join('\n', t) : null;
        return (p, targets);
    }

    private static Candidate JsonToCandidate(CandidateDto dto)
    {
        var ions = (dto.Fragments ?? Array.Empty<FragmentIonDto>())
            .Select(f => new FragmentIon(f.Mz, f.Intensity, f.Charge))
            .ToArray();
        return new Candidate
        {
            PrecursorId = dto.PrecursorId,
            StrippedSequence = dto.StrippedSequence,
            ModifiedSequence = dto.ModifiedSequence,
            PrecursorCharge = dto.PrecursorCharge,
            PrecursorMz = dto.PrecursorMz,
            RtStart = dto.RtStart,
            RtStop = dto.RtStop,
            RtApex = dto.RtApex,
            PrecursorQuantity = dto.PrecursorQuantity,
            PredictedPriority = dto.PredictedPriority,
            QValue = dto.QValue,
            ProteinQValue = dto.ProteinQValue,
            Proteotypic = dto.Proteotypic,
            ProteinGroup = dto.ProteinGroup,
            PeptideType = dto.PeptideType ?? "unique",
            Fragments = ions,
            Top4Fragments = Candidate.DeriveTopMz(ions, 4),
            Run = dto.Run ?? "cli",
        };
    }

    /// <summary>
    /// Overlay the DTO's set fields on <paramref name="defaults"/>; a null DTO
    /// gives the defaults. schedule / schedule-diann pass a bare
    /// <see cref="SchedulingParameters"/> (Core defaults); export-thermo passes
    /// <see cref="DiannThermoExport.AppDefaults"/>.
    /// </summary>
    internal static SchedulingParameters ToSchedulingParameters(SchedulingParametersDto? dto, SchedulingParameters defaults)
    {
        if (dto is null) return defaults;

        var d = defaults;
        return d with
        {
            Mode = dto.Mode ?? d.Mode,
            Objective = dto.Objective ?? d.Objective,
            IsolationWindowTh = dto.IsolationWindowTh ?? d.IsolationWindowTh,
            PrmIsolationWidthTh = dto.PrmIsolationWidthTh ?? d.PrmIsolationWidthTh,
            FragmentTolDa = dto.FragmentTolDa ?? d.FragmentTolDa,
            CycleBudget = dto.CycleBudget ?? d.CycleBudget,
            FiringPadSec = dto.FiringPadSec ?? d.FiringPadSec,
            RtBinMin = dto.RtBinMin ?? d.RtBinMin,
            QValueCutoff = dto.QValueCutoff ?? d.QValueCutoff,
            MinPeptidesPerProtein = dto.MinPeptidesPerProtein ?? d.MinPeptidesPerProtein,
            MaxPeptidesPerProtein = dto.MaxPeptidesPerProtein ?? d.MaxPeptidesPerProtein,
            PeptideRankCap = dto.PeptideRankCap ?? d.PeptideRankCap,
            MinTargetIonShare = dto.MinTargetIonShare ?? d.MinTargetIonShare,
            MaxTargetsPerSlot = dto.MaxTargetsPerSlot ?? d.MaxTargetsPerSlot,
            BackgroundFragmentClash = dto.BackgroundFragmentClash ?? d.BackgroundFragmentClash,
            ForcedPeptides = dto.ForcedPeptides is null
                ? d.ForcedPeptides
                : ForcedPeptides.Parse(string.Join('\n', dto.ForcedPeptides), out _),
            TargetMode = dto.TargetMode ?? d.TargetMode,
            ProteinRanking = dto.ProteinRanking ?? d.ProteinRanking,
            PeptideRanking = dto.PeptideRanking ?? d.PeptideRanking,
            ChargeHandling = dto.ChargeHandling ?? d.ChargeHandling,
            NormalizedCollisionEnergy = dto.NormalizedCollisionEnergy ?? d.NormalizedCollisionEnergy,
        };
    }

    private static ScheduleResponse ToResponse(ScheduleResult result)
    {
        return new ScheduleResponse
        {
            ScheduledIndices = result.ScheduledIndices,
            ScheduledSlotIds = result.ScheduledSlotIds,
            Slots = result.Slots.Select(s => new SlotDto
            {
                Id = s.Id,
                MzMin = s.MzMin,
                MzMax = s.MzMax,
                RtStart = s.RtStart,
                RtStop = s.RtStop,
                CoStart = s.CoStart,
                CoStop = s.CoStop,
                MemberIndices = s.MemberIndices.ToArray(),
            }).ToArray(),
            RtGrid = result.RtGrid,
            SlotCountCurve = result.SlotCountCurve,
            ProteinGroupsCovered = result.ProteinGroupsCovered,
        };
    }

    private static int Help()
    {
        PrintUsage();
        return 0;
    }

    private static int UnknownVerb(string verb)
    {
        Console.Error.WriteLine($"unknown verb '{verb}'");
        PrintUsage();
        return 64;
    }

    private static void PrintUsage()
    {
        Console.Error.WriteLine("Usage: SkylineCadenzaCli <verb> [args...]");
        Console.Error.WriteLine("Verbs:");
        Console.Error.WriteLine("  schedule INPUT.json OUTPUT.json   Run the scheduler.");
        Console.Error.WriteLine("                                     Pass '-' for stdin/stdout.");
        Console.Error.WriteLine("  schedule-diann REPORT.parquet PARAMS.json OUTPUT.json");
        Console.Error.WriteLine("                                     Ingest a DIA-NN report as the app does, then schedule.");
        Console.Error.WriteLine("  export-thermo REPORT.parquet PARAMS.json|- OUTPUT.csv");
        Console.Error.WriteLine("                                     Write the app's Thermo scheduled inclusion CSV.");
        Console.Error.WriteLine("                                     PARAMS '-' = the app's default settings.");
        Console.Error.WriteLine("                                     SUMMARY.json (optional) = what was scheduled.");
        Console.Error.WriteLine("                                     --skyline OUT.sky = also write the assay as a Skyline");
        Console.Error.WriteLine("                                     document (with its BLIB), via SkylineCmd");
        Console.Error.WriteLine("                                     (--skylinecmd PATH, else the newest installed).");
        Console.Error.WriteLine("  describe-diann REPORT.parquet PARAMS.json|- OUTPUT.json");
        Console.Error.WriteLine("                                     List the protein groups and precursors export-thermo");
        Console.Error.WriteLine("                                     would schedule from, by its names.");
        Console.Error.WriteLine("  help                              Show this help.");
    }
}

// ===== JSON DTOs =====================================================
// Kept separate from the Core types so the on-disk schema is stable
// across Core refactors. Property names are camelCased on the wire.

internal sealed class ScheduleRequest
{
    public CandidateDto[] Candidates { get; set; } = Array.Empty<CandidateDto>();
    public SchedulingParametersDto? Parameters { get; set; }
}

internal sealed class CandidateDto
{
    public string PrecursorId { get; set; } = string.Empty;
    public string StrippedSequence { get; set; } = string.Empty;
    public string ModifiedSequence { get; set; } = string.Empty;
    public int PrecursorCharge { get; set; }
    public double PrecursorMz { get; set; }
    public double RtStart { get; set; }
    public double RtStop { get; set; }
    public double RtApex { get; set; }
    public double PrecursorQuantity { get; set; }
    public double? PredictedPriority { get; set; }
    public double QValue { get; set; }
    public double ProteinQValue { get; set; } = double.NaN;
    public int Proteotypic { get; set; }
    public string ProteinGroup { get; set; } = string.Empty;
    public string? PeptideType { get; set; }
    public FragmentIonDto[]? Fragments { get; set; }
    public string? Run { get; set; }
}

internal sealed class FragmentIonDto
{
    public double Mz { get; set; }
    public double Intensity { get; set; }
    public int Charge { get; set; } = 1;
}

internal sealed class SchedulingParametersDto
{
    public AcquisitionMode? Mode { get; set; }
    public CoverageObjective? Objective { get; set; }
    public double? IsolationWindowTh { get; set; }
    public double? PrmIsolationWidthTh { get; set; }
    public double? FragmentTolDa { get; set; }
    public int? CycleBudget { get; set; }
    public double? FiringPadSec { get; set; }
    public double? RtBinMin { get; set; }
    public double? QValueCutoff { get; set; }
    public int? MinPeptidesPerProtein { get; set; }
    public int? MaxPeptidesPerProtein { get; set; }
    public int? PeptideRankCap { get; set; }
    public double? MinTargetIonShare { get; set; }
    public int? MaxTargetsPerSlot { get; set; }
    public bool? BackgroundFragmentClash { get; set; }
    /// <summary>Forced peptide entries, one per element, same syntax as the app's box.</summary>
    public string[]? ForcedPeptides { get; set; }
    /// <summary>
    /// export-thermo only: target proteins (accessions or gene names), one per
    /// element, same syntax as the app's box. Resolved against the report.
    /// </summary>
    public string[]? TargetProteins { get; set; }
    public TargetListMode? TargetMode { get; set; }
    public ProteinPriority? ProteinRanking { get; set; }
    public PeptidePriority? PeptideRanking { get; set; }
    public ChargeHandling? ChargeHandling { get; set; }
    public double? NormalizedCollisionEnergy { get; set; }
}

internal sealed class ScheduleResponse
{
    public int[] ScheduledIndices { get; set; } = Array.Empty<int>();
    public int[] ScheduledSlotIds { get; set; } = Array.Empty<int>();
    public SlotDto[] Slots { get; set; } = Array.Empty<SlotDto>();
    public double[] RtGrid { get; set; } = Array.Empty<double>();
    public int[] SlotCountCurve { get; set; } = Array.Empty<int>();
    public int ProteinGroupsCovered { get; set; }
}

internal sealed class SlotDto
{
    public int Id { get; set; }
    public double MzMin { get; set; }
    public double MzMax { get; set; }
    public double RtStart { get; set; }
    public double RtStop { get; set; }
    public double CoStart { get; set; }
    public double CoStop { get; set; }
    public int[] MemberIndices { get; set; } = Array.Empty<int>();
}

internal sealed class DiannScheduleResponse
{
    public int PoolPrecursors { get; set; }
    public int ProteinGroupsCovered { get; set; }
    public int Slots { get; set; }
    public ScheduledTargetDto[] Targets { get; set; } = Array.Empty<ScheduledTargetDto>();
    public ForcedOutcomeDto[] Forced { get; set; } = Array.Empty<ForcedOutcomeDto>();
}

internal sealed class ReportDescription
{
    public double QValueCutoff { get; set; }
    public int ProteinGroups { get; set; }
    public int Precursors { get; set; }
    public IReadOnlyList<ReportCatalog.ProteinEntry> Proteins { get; set; } = Array.Empty<ReportCatalog.ProteinEntry>();
}

internal sealed class ForcedOutcomeDto
{
    public string Entry { get; set; } = string.Empty;
    public string Status { get; set; } = string.Empty;
    public string? PrecursorId { get; set; }
    public bool UsedFallbackCharge { get; set; }
}

internal sealed class ScheduledTargetDto
{
    public string PrecursorId { get; set; } = string.Empty;
    public string StrippedSequence { get; set; } = string.Empty;
    public int PrecursorCharge { get; set; }
    public string ProteinGroup { get; set; } = string.Empty;
    public double PrecursorMz { get; set; }
    public double PrecursorQuantity { get; set; }
    public int SlotId { get; set; }
    public double SlotMzMin { get; set; }
    public double SlotMzMax { get; set; }
    public double SlotRtStart { get; set; }
    public double SlotRtStop { get; set; }
    public int SlotSize { get; set; }
}
