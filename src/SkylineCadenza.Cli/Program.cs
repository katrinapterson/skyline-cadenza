using System.Text.Json;
using System.Text.Json.Serialization;
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
            ProteinRanking = dto.ProteinRanking ?? ProteinPriority.SummedIntensity,
            PeptideRanking = dto.PeptideRanking ?? PeptidePriority.PrecursorQuantity,
            ChargeHandling = dto.ChargeHandling ?? ChargeHandling.SameChargePerSlot,
            NormalizedCollisionEnergy = dto.NormalizedCollisionEnergy ?? 28.0,
            Packing = dto.Packing ?? WindowPacking.Legacy,
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
    public ProteinPriority? ProteinRanking { get; set; }
    public PeptidePriority? PeptideRanking { get; set; }
    public ChargeHandling? ChargeHandling { get; set; }
    public double? NormalizedCollisionEnergy { get; set; }
    public WindowPacking? Packing { get; set; }
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
