using System.Globalization;
using Parquet;
using Parquet.Data;
using Parquet.Schema;
using SkylineCadenza.Cli;
using SkylineCadenza.Core.Ingest;
using SkylineCadenza.Core.Output;
using SkylineCadenza.Core.Parsimony;
using SkylineCadenza.Core.Scheduling;
using SkylineCadenza.Core.SkylineRpc;
using Xunit;
using Xunit.Abstractions;

namespace SkylineCadenza.Tests;

/// <summary>
/// <c>SkylineCadenzaCli export-thermo</c> must write the same CSV the app's
/// "Export Thermo CSV..." writes after an Empirical Load with the same
/// settings. The reference here composes the Core calls the way
/// MainViewModel.LoadFromDiannAsync + ExportThermoCsv do.
/// </summary>
public class ExportThermoCliTests : IDisposable
{
    // A real Stellar GPF DIA-NN report; the comparison on it is skipped when absent.
    private const string RealReportPath =
        @"C:\Users\MacCoss Lab\Katrina\cadenza\stellar_08262026\GPF_astral-IDs_carafe\diann_project\report.parquet";

    private readonly ITestOutputHelper _output;
    private readonly string _dir;

    public ExportThermoCliTests(ITestOutputHelper output)
    {
        _output = output;
        _dir = Path.Combine(Path.GetTempPath(), "cadenza-export-thermo-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    public void Dispose()
    {
        try { Directory.Delete(_dir, recursive: true); } catch { /* best effort */ }
    }

    [Fact]
    public void App_defaults_differ_from_core_only_where_the_app_does()
    {
        var core = new SchedulingParameters();
        var app = DiannThermoExport.AppDefaults;
        // MainViewModel: _minTargetIonShare = 0.015, BackgroundFragmentClash = true always.
        Assert.Equal(0.015, app.MinTargetIonShare);
        Assert.True(app.BackgroundFragmentClash);
        // (TargetProteins is a fresh empty set per instance; share it so record equality compares values.)
        Assert.Equal(core with { MinTargetIonShare = 0.015, BackgroundFragmentClash = true, TargetProteins = app.TargetProteins }, app);
    }

    [Fact]
    public void Dash_or_empty_params_give_the_app_defaults()
    {
        Assert.Equal(DiannThermoExport.AppDefaults, Program.ReadExportParameters(null).Parameters);
        Assert.Equal(DiannThermoExport.AppDefaults, Program.ReadExportParameters("{}").Parameters);
        Assert.Null(Program.ReadExportParameters("{}").TargetProteinsText);
    }

    [Fact]
    public void Params_json_overrides_only_the_fields_it_sets()
    {
        var (p, targets) = Program.ReadExportParameters("""
            {
              "mode": "Prm",
              "cycleBudget": 60,
              "normalizedCollisionEnergy": 30,
              "backgroundFragmentClash": false,
              "targetProteins": ["P00001", "GENE2"],
              "targetMode": "Exclusive"
            }
            """);
        var expected = DiannThermoExport.AppDefaults with
        {
            Mode = AcquisitionMode.Prm,
            CycleBudget = 60,
            NormalizedCollisionEnergy = 30,
            BackgroundFragmentClash = false,
            TargetMode = TargetListMode.Exclusive,
        };
        Assert.Equal(expected, p);
        Assert.Equal("P00001\nGENE2", targets);
    }

    [Fact]
    public void Other_verbs_keep_core_defaults()
    {
        var core = new SchedulingParameters();
        Assert.Same(core, Program.ToSchedulingParameters(null, core));
        Assert.Equal(core, Program.ToSchedulingParameters(new SchedulingParametersDto(), core));
        Assert.Equal(0.0, Program.ToSchedulingParameters(new SchedulingParametersDto(), new SchedulingParameters()).MinTargetIonShare);
    }

    [Fact]
    public async Task Missing_report_fails_with_a_clear_exit_code()
    {
        int rc = await Program.RunExportThermoAsync(new[]
        {
            Path.Combine(_dir, "nope.parquet"), "-", Path.Combine(_dir, "out.csv"),
        });
        Assert.Equal(66, rc);
        Assert.False(File.Exists(Path.Combine(_dir, "out.csv")));
        Assert.Equal(64, await Program.RunExportThermoAsync(new[] { "only-one-arg" }));
    }

    [Fact]
    public async Task Synthetic_report_matches_the_app_export_in_mtm()
    {
        string report = Path.Combine(_dir, "report.parquet");
        await WriteSyntheticReportAsync(report);
        string csvPath = Path.Combine(_dir, "sub", "cadenza_targets.csv");

        int rc = await Program.RunExportThermoAsync(new[] { report, "-", csvPath });
        Assert.Equal(0, rc);

        var (expected, schedule) = await AppExportAsync(report, DiannThermoExport.AppDefaults);
        string actual = File.ReadAllText(csvPath);
        Assert.Equal(expected, actual);

        var lines = actual.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        Assert.StartsWith("Compound,Formula,Adduct,m/z,z,t start (min),t stop (min),Isolation Window (m/z),HCD Collision Energy", lines[0]);
        // MTM: one row per non-empty slot.
        Assert.Equal(schedule.Slots.Count(s => s.MemberIndices.Count > 0), lines.Length - 1);
        Assert.True(lines.Length - 1 > 0);
        // Multiplexed slots exist in this fixture and are written as one joined row.
        Assert.Contains(lines, l => l.Contains(" | "));
    }

    [Fact]
    public async Task Synthetic_report_matches_the_app_export_in_prm_with_params_file()
    {
        string report = Path.Combine(_dir, "report.parquet");
        await WriteSyntheticReportAsync(report);
        string paramsPath = Path.Combine(_dir, "params.json");
        File.WriteAllText(paramsPath, """{ "mode": "Prm", "normalizedCollisionEnergy": 27.5 }""");
        string csvPath = Path.Combine(_dir, "prm.csv");

        Assert.Equal(0, await Program.RunExportThermoAsync(new[] { report, paramsPath, csvPath }));

        var p = DiannThermoExport.AppDefaults with { Mode = AcquisitionMode.Prm, NormalizedCollisionEnergy = 27.5 };
        var (expected, schedule) = await AppExportAsync(report, p);
        string actual = File.ReadAllText(csvPath);
        Assert.Equal(expected, actual);
        var rows = actual.Split('\n', StringSplitOptions.RemoveEmptyEntries).Skip(1).ToList();
        Assert.Equal(schedule.ScheduledIndices.Length, rows.Count);
        Assert.All(rows, r => Assert.EndsWith("0.7000,27.5", r.TrimEnd('\r')));
    }

    [Fact]
    public async Task Real_report_matches_the_app_export()
    {
        if (!File.Exists(RealReportPath))
        {
            _output.WriteLine($"Reference report missing - skipping. {RealReportPath}");
            return;
        }
        string csvPath = Path.Combine(_dir, "real.csv");
        Assert.Equal(0, await Program.RunExportThermoAsync(new[] { RealReportPath, "-", csvPath }));
        var (expected, schedule) = await AppExportAsync(RealReportPath, DiannThermoExport.AppDefaults);
        Assert.Equal(expected, File.ReadAllText(csvPath));
        _output.WriteLine($"scheduled {schedule.ScheduledIndices.Length:n0} in {schedule.Slots.Length:n0} slots");
    }

    /// <summary>
    /// The app's path: Empirical Load (no Carafe library, no target list)
    /// then "Export Thermo CSV..." on the resulting schedule.
    /// </summary>
    private static async Task<(string Csv, ScheduleResult Schedule)> AppExportAsync(string report, SchedulingParameters p)
    {
        var rows = await DiannParquetReader.LoadAsync(report, p.QValueCutoff);
        var parsimony = ParsimonyEngine.Assign(DiannParquetReader.BuildPeptideProteinMap(rows));
        var candidates = CandidateBuilder.Build(rows, parsimony, new Dictionary<FragmentKey, FragmentIon[]>());
        var schedule = Scheduler.Run(candidates, p);
        return (ThermoCsvWriter.Build(candidates, schedule, p), schedule);
    }

    /// <summary>
    /// A small DIA-NN-shaped report: 4 proteins x 4 charge-2 peptides with
    /// overlapping RT windows and close m/z, so MTM multiplexes some slots.
    /// Fragments are distinct per precursor so no fragment clash fires.
    /// </summary>
    private static async Task WriteSyntheticReportAsync(string path)
    {
        var run = new List<string>();
        var precursorId = new List<string>();
        var modSeq = new List<string>();
        var stripSeq = new List<string>();
        var charge = new List<long>();
        var mz = new List<double>();
        var rt = new List<double>();
        var rtStart = new List<double>();
        var rtStop = new List<double>();
        var quantity = new List<double>();
        var proteotypic = new List<long>();
        var proteinIds = new List<string>();
        var genes = new List<string>();
        var qValue = new List<double>();
        var decoy = new List<long>();
        var pgQ = new List<double>();
        var fragIds = Enumerable.Range(0, 4).Select(_ => new List<string>()).ToArray();
        var fragQty = Enumerable.Range(0, 4).Select(_ => new List<double>()).ToArray();

        const string aa = "ACDEFGHIKLMNPQRSTVWY";
        var inv = CultureInfo.InvariantCulture;
        int n = 0;
        for (int prot = 1; prot <= 4; prot++)
        {
            for (int pep = 0; pep < 4; pep++, n++)
            {
                // Unique tryptic-looking sequence per precursor.
                string seq = $"{aa[n % 20]}{aa[(n * 7 + 3) % 20]}{aa[(n * 3 + 1) % 20]}LVPEGTK";
                run.Add("run1");
                precursorId.Add(seq + "2");
                modSeq.Add(seq);
                stripSeq.Add(seq);
                charge.Add(2);
                // Pairs within ~1 Th of each other so MTM can share a window.
                double precursorMz = 500.0 + (n / 2) * 25.0 + (n % 2) * 0.9;
                mz.Add(precursorMz);
                double apex = 10.0 + (n / 2 % 5) * 0.4; // each m/z pair co-elutes
                rt.Add(apex);
                rtStart.Add(apex - 0.15);
                rtStop.Add(apex + 0.15);
                quantity.Add(1.0e6 * (1 + (n * 13) % 7));
                proteotypic.Add(1);
                proteinIds.Add($"P0000{prot}");
                genes.Add($"GENE{prot}");
                qValue.Add(0.001);
                decoy.Add(0);
                pgQ.Add(0.001);
                for (int k = 0; k < 4; k++)
                {
                    double fmz = 300.0 + n * 37.0 + k * 101.3;
                    fragIds[k].Add($"y{k + 3}^1/{fmz.ToString("0.000000", inv)}");
                    fragQty[k].Add(1.0e5 * (4 - k));
                }
            }
        }
        // A decoy and a high q-value row, both of which the reader must drop.
        run.Add("run1"); precursorId.Add("DECOYK2"); modSeq.Add("DECOYK"); stripSeq.Add("DECOYK");
        charge.Add(2); mz.Add(600); rt.Add(10); rtStart.Add(9.9); rtStop.Add(10.1); quantity.Add(1e9);
        proteotypic.Add(1); proteinIds.Add("P00001"); genes.Add("GENE1"); qValue.Add(0.001); decoy.Add(1); pgQ.Add(0.001);
        run.Add("run1"); precursorId.Add("BADQK2"); modSeq.Add("BADQK"); stripSeq.Add("BADQK");
        charge.Add(2); mz.Add(610); rt.Add(10); rtStart.Add(9.9); rtStop.Add(10.1); quantity.Add(1e9);
        proteotypic.Add(1); proteinIds.Add("P00002"); genes.Add("GENE2"); qValue.Add(0.5); decoy.Add(0); pgQ.Add(0.001);
        for (int k = 0; k < 4; k++)
        {
            fragIds[k].Add($"y{k + 3}^1/{(900.0 + k).ToString("0.000000", inv)}"); fragQty[k].Add(1e5);
            fragIds[k].Add($"y{k + 3}^1/{(950.0 + k).ToString("0.000000", inv)}"); fragQty[k].Add(1e5);
        }

        var columns = new List<(DataField Field, Array Data)>
        {
            (new DataField<string>("Run"), run.ToArray()),
            (new DataField<string>("Precursor.Id"), precursorId.ToArray()),
            (new DataField<string>("Modified.Sequence"), modSeq.ToArray()),
            (new DataField<string>("Stripped.Sequence"), stripSeq.ToArray()),
            (new DataField<long>("Precursor.Charge"), charge.ToArray()),
            (new DataField<double>("Precursor.Mz"), mz.ToArray()),
            (new DataField<double>("RT"), rt.ToArray()),
            (new DataField<double>("RT.Start"), rtStart.ToArray()),
            (new DataField<double>("RT.Stop"), rtStop.ToArray()),
            (new DataField<double>("Precursor.Quantity"), quantity.ToArray()),
            (new DataField<long>("Proteotypic"), proteotypic.ToArray()),
            (new DataField<string>("Protein.Ids"), proteinIds.ToArray()),
            (new DataField<string>("Genes"), genes.ToArray()),
            (new DataField<double>("Q.Value"), qValue.ToArray()),
            (new DataField<long>("Decoy"), decoy.ToArray()),
            (new DataField<double>("PG.Q.Value"), pgQ.ToArray()),
        };
        for (int k = 0; k < 4; k++)
        {
            columns.Add((new DataField<string>($"Fr.{k}.Id"), fragIds[k].ToArray()));
            columns.Add((new DataField<double>($"Fr.{k}.Quantity"), fragQty[k].ToArray()));
        }

        var schema = new ParquetSchema(columns.Select(c => (Field)c.Field).ToArray());
        using var stream = File.Create(path);
        using var writer = await ParquetWriter.CreateAsync(schema, stream);
        using var rg = writer.CreateRowGroup();
        foreach (var (field, data) in columns)
            await rg.WriteColumnAsync(new DataColumn(field, data));
    }
}
