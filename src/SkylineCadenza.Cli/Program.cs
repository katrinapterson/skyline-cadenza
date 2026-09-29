using System.Text.Json;
using System.Text.Json.Serialization;
using SkylineCadenza.Core.Ingest;
using SkylineCadenza.Core.Parsimony;
using SkylineCadenza.Core.Scheduling;

namespace SkylineCadenza.Cli;

/// <summary>
/// Headless companion to the WPF app. Reads a candidate pool + a
/// <see cref="SchedulingParameters"/> from JSON, runs
/// <see cref="Scheduler.Run"/>, and writes the schedule back out as
/// JSON. Used by the algorithm-comparison Jupyter notebook so that the
/// realistic plots reflect Cadenza's real C# scheduler bit-for-bit
/// rather than a Python re-implementation.
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
        var parameters = ToSchedulingParameters(request.Parameters);
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
        var parameters = ToSchedulingParameters(paramDto);
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

    private static SchedulingParameters ToSchedulingParameters(SchedulingParametersDto? dto)
    {
        // Empty / null DTO -> stock defaults.
        if (dto is null) return new SchedulingParameters();

        var p = new SchedulingParameters
        {
            Mode = dto.Mode ?? AcquisitionMode.Mtm,
            Objective = dto.Objective ?? CoverageObjective.MaximizeProteins,
            IsolationWindowTh = dto.IsolationWindowTh ?? 3.0,
            PrmIsolationWidthTh = dto.PrmIsolationWidthTh ?? 0.7,
            FragmentTolDa = dto.FragmentTolDa ?? 0.5,
            CycleBudget = dto.CycleBudget ?? 100,
            FiringPadSec = dto.FiringPadSec ?? 15.0,
            RtBinMin = dto.RtBinMin ?? 0.05,
            QValueCutoff = dto.QValueCutoff ?? 0.01,
            MinPeptidesPerProtein = dto.MinPeptidesPerProtein ?? 1,
            MaxPeptidesPerProtein = dto.MaxPeptidesPerProtein ?? 5,
            PeptideRankCap = dto.PeptideRankCap ?? 0,
            MinTargetIonShare = dto.MinTargetIonShare ?? 0.0,
            MaxTargetsPerSlot = dto.MaxTargetsPerSlot ?? 0,
            BackgroundFragmentClash = dto.BackgroundFragmentClash ?? false,
            ForcedPeptides = ForcedPeptides.Parse(string.Join('\n', dto.ForcedPeptides ?? Array.Empty<string>()), out _),
            ProteinRanking = dto.ProteinRanking ?? ProteinPriority.SummedIntensity,
            PeptideRanking = dto.PeptideRanking ?? PeptidePriority.PrecursorQuantity,
            ChargeHandling = dto.ChargeHandling ?? ChargeHandling.SameChargePerSlot,
            NormalizedCollisionEnergy = dto.NormalizedCollisionEnergy ?? 28.0,
        };
        return p;
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
