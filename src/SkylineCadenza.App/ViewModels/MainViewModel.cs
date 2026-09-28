using System.IO;
using System.Threading;
using System.Threading.Tasks;
using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;
using Microsoft.Win32;
using SkylineCadenza.Core.Ingest;
using SkylineCadenza.Core.Output;
using SkylineCadenza.Core.Parsimony;
using SkylineCadenza.Core.Scheduling;
using SkylineCadenza.Core.SkylineRpc;
using CoverageCurves = SkylineCadenza.Core.Scheduling.CoverageCurves;

namespace SkylineCadenza.App.ViewModels;

/// <summary>
/// Root view model: input file paths, scheduling parameters (live-bound),
/// async commands to ingest + schedule, and an observable
/// <see cref="ScheduleResult"/> that the plot views render from.
/// </summary>
public partial class MainViewModel : ObservableObject
{
    [ObservableProperty] private string? _reportPath;
    [ObservableProperty] private string? _carafeTsvPath;

    /// <summary>
    /// Optional Bromo peptide-ranking overlay for Carafe-only mode: a single
    /// Bromo <c>peptide_rankings.tsv</c>. When set, peptide priority
    /// auto-switches to <see cref="PeptidePriority.PredictedRanking"/> at load.
    /// (Carafe detectability predictions go in <see cref="DetectabilityPath"/>.)
    /// </summary>
    [ObservableProperty] private string? _predictedRankingPath;

    /// <summary>
    /// Synthesised peak half-width (seconds) used on the Carafe path, where the
    /// predicted library carries only an apex. The candidate's unpadded peak is
    /// <c>apex +/- this</c>; the scheduler then adds <see cref="FiringPadSec"/>
    /// on top. Ignored on the DIA-NN path (real RT.Start/RT.Stop are used).
    /// 6 s = the legacy 0.10 min default.
    /// </summary>
    [ObservableProperty] private double _peakHalfWidthSec = 6.0;

    /// <summary>Convenience accessor: <see cref="PeakHalfWidthSec"/> in minutes.</summary>
    public double PeakHalfWidthMin => PeakHalfWidthSec / 60.0;

    /// <summary>
    /// Optional Carafe detectability prediction files (one or more
    /// <c>*_detectability_pred.tsv</c>, or a folder) used as a candidate FILTER
    /// on the Carafe path - separate from the ranking overlay above, so the user
    /// can rank with Bromo and filter with Carafe detectability at the same time.
    /// </summary>
    [ObservableProperty] private string? _detectabilityPath;

    /// <summary>
    /// Minimum Carafe detectability a precursor must have to be kept (Carafe
    /// path). Applied only when <see cref="DetectabilityPath"/> is set;
    /// precursors with no prediction are kept. Default 0.5 (empirically drops
    /// peptides ~90% unlikely to be observed while retaining most detections).
    /// </summary>
    [ObservableProperty] private double _detectabilityCutoff = 0.5;

    /// <summary>
    /// When a DIA-NN report and a Carafe library are both loaded, schedule the
    /// confirmed DIA-NN targets first, then backfill leftover schedule capacity
    /// with Carafe-predicted targets (deduped against the confirmed ones on
    /// peptide + charge). Off = the Carafe library only overrides DIA-NN
    /// fragments and adds no precursors. Takes effect on the next Load.
    /// </summary>
    [ObservableProperty] private bool _backfillWithPredictions;

    // Acquisition + load balancing
    [ObservableProperty]
    [NotifyPropertyChangedFor(nameof(SlotsLabel))]
    private AcquisitionMode _mode = AcquisitionMode.Mtm;

    /// <summary>Header text for the slots metric: "MTM slots" or "PRM slots".</summary>
    public string SlotsLabel => Mode == AcquisitionMode.Prm ? "PRM slots" : "MTM slots";
    [ObservableProperty] private CoverageObjective _objective = CoverageObjective.MaximizeProteins;

    // MTM constraints
    [ObservableProperty] private double _isolationWindowTh = 3.0;
    [ObservableProperty] private double _fragmentTolDa = 0.5;

    /// <summary>
    /// Minimum share of the co-isolated ion current that EVERY target in a
    /// slot must hold while it elutes. 0 disables the check. Defaults to 1.5%,
    /// the break-even floor on the acquired 08/03 run once truncated-FASTA
    /// artefact targets are excluded from the label: it removes 0.98 targets
    /// DIA-NN never re-identified per one it does cost, against 0.82 at 2%.
    /// (2% looked break-even before the artefacts were removed, because the
    /// gate was getting credit for excluding targets that were never real
    /// peptides.) Costs 1.7% of targets and 1.8% of driver protein groups.
    /// </summary>
    [ObservableProperty] private double _minTargetIonShare = 0.015;

    /// <summary>
    /// Include passengers when pushing to Skyline. Passengers already ride in
    /// the scheduled windows, so this changes only what the document extracts,
    /// never the acquisition: the exported Thermo method is identical either
    /// way. Off by default so the pushed document holds exactly what was
    /// targeted.
    /// </summary>
    [ObservableProperty] private bool _includePassengers;

    // PRM quad isolation width (also the minimum slot width for solo MTM slots).
    [ObservableProperty] private double _prmIsolationWidthTh = 0.7;

    // How isolation windows may grow relative to each other. Legacy keeps the
    // historical behaviour (and the golden fixtures); Disjoint stops a window
    // expanding across a clashing precursor and forbids concurrent windows from
    // overlapping in m/z. See docs/scheduler-window-scoping.md.

    // Per-slot charge constraint + Thermo CSV NCE.
    [ObservableProperty] private ChargeHandling _chargeHandling = ChargeHandling.SameChargePerSlot;
    [ObservableProperty] private double _normalizedCollisionEnergy = 28.0;

    // Cycle budget
    [ObservableProperty] private int _cycleBudget = 100;
    [ObservableProperty] private double _firingPadSec = 15.0;

    // Library filter
    [ObservableProperty] private double _qValueCutoff = 0.01;

    // Per-protein bounds + priority knobs
    [ObservableProperty] private int _minPeptidesPerProtein = 1;
    [ObservableProperty] private int _maxPeptidesPerProtein = 5;

    /// <summary>
    /// Only the top N precursors of each protein group (by the active
    /// <see cref="PeptideRanking"/>, per source tier) are eligible for
    /// scheduling. 0 = no cap.
    /// </summary>
    [ObservableProperty] private int _peptideRankCap;
    [ObservableProperty] private ProteinPriority _proteinRanking = ProteinPriority.SummedIntensity;
    [ObservableProperty] private PeptidePriority _peptideRanking = PeptidePriority.PrecursorQuantity;

    // Target protein list
    [ObservableProperty] private string _targetProteinsText = string.Empty;
    [ObservableProperty] private TargetListMode _targetMode = TargetListMode.FirstThenFill;
    [ObservableProperty] private int _targetProteinsParsedCount;

    // Status + result
    [ObservableProperty] private string _statusMessage = "Open a DIA-NN report.parquet to begin.";
    [ObservableProperty] private bool _isBusy;

    // Library-level stats (everything that passed the q-value filter in the
    // report, before scheduling).
    [ObservableProperty] private int _libraryPrecursors;
    [ObservableProperty] private int _libraryPeptides;
    [ObservableProperty] private int _libraryProteinGroups;
    [ObservableProperty] private int _libraryCarafeMatched;

    // Heatmap view toggle.
    [ObservableProperty]
    [NotifyPropertyChangedFor(nameof(HeatmapTitle))]
    private HeatmapView _heatmapView = HeatmapView.FullLibrary;

    public string HeatmapTitle => HeatmapView == HeatmapView.FullLibrary
        ? "Precursors per spectrum (full library)"
        : "Precursors per spectrum (scheduled assay)";

    [ObservableProperty] private ScheduleResult? _scheduleResult;
    [ObservableProperty] private int _scheduledPrecursors;
    [ObservableProperty] private int _scheduledPeptides;
    [ObservableProperty] private int _proteinsCovered;
    [ObservableProperty] private int _slotsTotal;
    [ObservableProperty] private int _slotsShared;
    [ObservableProperty] private int _peakLoad;
    [ObservableProperty] private double _medianPeptidesPerProtein;

    /// <summary>Per-protein summary used by the protein-coverage plot.</summary>
    public List<ProteinCoverageRow> ProteinCoverage { get; } = new();

    /// <summary>Per-run coverage curves (one per distinct DIA-NN run name).</summary>
    public CoverageCurves.CurveSeries[] PerRunCurves { get; private set; } = Array.Empty<CoverageCurves.CurveSeries>();

    /// <summary>m/z x RT heatmap of the full library.</summary>
    public CoverageCurves.Heatmap? FullLibraryHeatmap { get; private set; }

    /// <summary>m/z x RT heatmap of just the scheduled precursors.</summary>
    public CoverageCurves.Heatmap? ScheduledHeatmap { get; private set; }

    /// <summary>The heatmap currently selected by <see cref="HeatmapView"/>.</summary>
    public CoverageCurves.Heatmap? Heatmap =>
        HeatmapView == HeatmapView.FullLibrary ? FullLibraryHeatmap : ScheduledHeatmap;

    /// <summary>Raised after a successful reschedule so views refresh.</summary>
    public event Action? ScheduleUpdated;

    /// <summary>Raised after a fresh ingest. Library-level plots refresh from this.</summary>
    public event Action? LibraryLoaded;

    private List<Candidate>? _candidates;
    private CancellationTokenSource? _rescheduleCts;
    private Dictionary<string, IReadOnlyList<string>>? _geneToAccession;
    private HashSet<string>? _knownAccessions;
    private SkylineSession? _skylineSession;

    // Cached Carafe-path inputs so the peak-half-width slider can rebuild
    // candidates (RtStart/RtStop are baked in at build time) and reschedule
    // live, without re-streaming the Carafe library. Null unless the active
    // pool came from the Carafe workflow.
    private List<CarafePrecursor>? _carafePrecursors;
    private IReadOnlyDictionary<string, (string Group, string Type)>? _carafeParsimony;
    private PredictedRankingReader.Result? _carafePredicted;
    private PredictedRankingReader.Result? _carafeDetectability;
    // Combined DIA-NN + Carafe backfill: when non-null, the cached Carafe
    // precursors above are the PREDICTED tier and these are the EXPERIMENTAL
    // tier's inputs, so the sliders rebuild the whole two-tier pool live.
    private List<DiannRow>? _combinedDiannRows;
    private IReadOnlyDictionary<FragmentKey, FragmentIon[]>? _combinedFrags;
    private System.Threading.Timer? _skylinePoll;
    private int _skylinePollFailures;
    private CancellationTokenSource? _skylineLoadCts;

    /// <summary>True while a long Skyline operation is interruptible via <see cref="CancelSkylineLoadCommand"/>.</summary>
    [ObservableProperty] private bool _skylineLoadCancellable;

    /// <summary>
    /// Optional override for the Skyline RT column. When set, the loader
    /// skips its auto-probe and uses this column name verbatim. Useful
    /// for documents whose RT column isn't in the loader's default list.
    /// </summary>
    [ObservableProperty] private string _skylineRtColumnOverride = string.Empty;

    [RelayCommand]
    private void CancelSkylineLoad()
    {
        _skylineLoadCts?.Cancel();
        StatusMessage = "Skyline load cancellation requested...";
        SkylineStatus = "Skyline: load cancellation requested.";
    }

    // Skyline connection state (bound to the toolbar in the UI).
    [ObservableProperty] private bool _skylineConnected;
    [ObservableProperty] private string _skylineStatus = "Skyline: not connected.";
    [ObservableProperty] private string? _skylineDocumentPath;

    public MainViewModel()
    {
        // Try to attach to a running Skyline instance on launch. Done on a
        // background thread because NamedPipeClientStream.Connect() blocks
        // for up to 5 s during discovery; running it inline would freeze
        // the splash.
        _ = Task.Run(TryConnectSkyline);
    }

    private void TryConnectSkyline()
    {
        try
        {
            _skylineSession = SkylineSession.FromArguments(App.LaunchArgs);
            // Don't open the pipe yet - Skyline closes it after each
            // request, so we connect lazily inside each Execute() call.
            // Kick off a 2 s polling loop to track document changes.
            _skylinePoll = new System.Threading.Timer(_ => PollSkyline(),
                null, TimeSpan.Zero, TimeSpan.FromSeconds(2));
        }
        catch (Exception ex)
        {
            System.Windows.Application.Current?.Dispatcher.Invoke(() =>
            {
                SkylineConnected = false;
                SkylineStatus = $"Skyline: not connected ({ex.GetType().Name}: {ex.Message}).";
            });
        }
    }

    private void PollSkyline()
    {
        if (_skylineSession is null) return;
        try
        {
            // Connect-per-call: each Skyline RPC opens a fresh pipe and
            // disposes after the response. Reusing a single pipe across
            // calls fails because Skyline closes its side between requests.
            var st = _skylineSession.Execute(c => c.GetDocumentStatus());
            _skylinePollFailures = 0;
            System.Windows.Application.Current?.Dispatcher.Invoke(() =>
            {
                SkylineConnected = true;
                SkylineDocumentPath = st?.DocumentPath;
                SkylineStatus = st is null
                    ? "Skyline: connected, no document open."
                    : $"Skyline: {st.DocumentPath ?? "(unsaved)"} - {st.Molecules:n0} peptides, {st.Precursors:n0} precursors.";
            });
        }
        catch (System.Text.Json.JsonException jex)
        {
            // Most often this means Skyline sent back an empty or
            // truncated message - it happens transiently when the pipe is
            // mid-handshake. Surface the exact reader message so we can
            // tell apart "empty response" from "unexpected token".
            _skylinePollFailures++;
            if (_skylinePollFailures >= 3)
            {
                System.Windows.Application.Current?.Dispatcher.Invoke(() =>
                {
                    SkylineConnected = false;
                    SkylineStatus = $"Skyline: response not parseable as JSON ({jex.Message}).";
                });
            }
        }
        catch (Exception ex)
        {
            // Skyline closes pipes between requests, so a transient
            // IOException right after a successful poll is normal. Only
            // flip to "connection lost" after a couple of failures.
            _skylinePollFailures++;
            if (_skylinePollFailures >= 3)
            {
                System.Windows.Application.Current?.Dispatcher.Invoke(() =>
                {
                    SkylineConnected = false;
                    SkylineStatus = $"Skyline: connection lost ({ex.GetType().Name}: {ex.Message}).";
                });
            }
        }
    }

    [RelayCommand]
    private void BrowseReport()
    {
        var dlg = new OpenFileDialog
        {
            Filter = "DIA-NN report (*.parquet)|*.parquet|All files (*.*)|*.*",
            Title = "Select DIA-NN report.parquet",
        };
        if (dlg.ShowDialog() == true)
            ReportPath = dlg.FileName;
    }

    [RelayCommand]
    private void BrowseCarafe()
    {
        var dlg = new OpenFileDialog
        {
            Filter = "Carafe library (*.tsv)|*.tsv|All files (*.*)|*.*",
            Title = "Select Carafe spectral library TSV",
        };
        if (dlg.ShowDialog() == true)
            CarafeTsvPath = dlg.FileName;
    }

    [RelayCommand]
    private void BrowsePredictedRanking()
    {
        var dlg = new OpenFileDialog
        {
            Filter = "Bromo peptide ranking (*.tsv)|*.tsv|All files (*.*)|*.*",
            Title = "Select a Bromo peptide_rankings.tsv",
        };
        if (dlg.ShowDialog() == true)
            PredictedRankingPath = dlg.FileName;
    }

    [RelayCommand]
    private void BrowseDetectability()
    {
        var dlg = new OpenFileDialog
        {
            Filter = "Carafe detectability (*.tsv)|*.tsv|All files (*.*)|*.*",
            Title = "Select one or more Carafe *_detectability_pred.tsv files",
            Multiselect = true,
        };
        if (dlg.ShowDialog() == true)
            DetectabilityPath = string.Join(Environment.NewLine, dlg.FileNames);
    }

    /// <summary>
    /// Split the predicted-ranking box into individual paths. Each non-blank
    /// line is one file or folder; supports a single path, a folder, or several
    /// multi-selected detectability shards.
    /// </summary>
    private static List<string> ParsePredictedRankingPaths(string? text) =>
        string.IsNullOrWhiteSpace(text)
            ? new List<string>()
            : text.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries
                    | StringSplitOptions.TrimEntries)
                .ToList();

    [RelayCommand]
    private void ExportThermoCsv()
    {
        if (_candidates is null || ScheduleResult is null)
        {
            StatusMessage = "Schedule something first - no targets to export.";
            return;
        }
        var dlg = new Microsoft.Win32.SaveFileDialog
        {
            Filter = "Thermo inclusion list (*.csv)|*.csv|All files (*.*)|*.*",
            FileName = "cadenza_targets.csv",
            Title = "Export Thermo scheduled inclusion CSV",
        };
        if (dlg.ShowDialog() != true) return;
        try
        {
            var csv = ThermoCsvWriter.Build(_candidates, ScheduleResult, CurrentParameters());
            File.WriteAllText(dlg.FileName, csv);
            StatusMessage = $"Wrote {dlg.FileName}. Import this in the Thermo Method Editor (Import Scheduled List).";
        }
        catch (Exception ex)
        {
            StatusMessage = $"Failed to write CSV: {ex.Message}";
        }
    }

    [RelayCommand]
    private async Task UpdateSkylineAsync()
    {
        if (_candidates is null || ScheduleResult is null)
        {
            StatusMessage = "Schedule something first - no targets to push.";
            return;
        }
        if (_skylineSession is null)
        {
            StatusMessage = "Skyline: not connected. Launch this tool from Skyline's Tools menu.";
            return;
        }

        // Two-step push:
        //   1. Write a self-contained BiblioSpec library (one reference
        //      spectrum per peptide with top-6 library fragments + the
        //      per-replicate RT firing window) and register it as a
        //      peptide-settings library in the running document. Skyline
        //      uses this BLIB for chromatogram extraction and as the
        //      canonical RT / fragment-intensity record of the assay.
        //
        //   2. Push a peptide-style transition list so the document's
        //      target tree is populated immediately. Fragment Product
        //      Charge values come from FragmentIon.Charge directly, so
        //      the +2 fragments that triggered "no matching product ion"
        //      errors in the previous all-+1-hardcoded build are now
        //      reported correctly.
        string? blibPath = null;
        string? tlistPath = null;
        try
        {
            // Local time so the timestamp on the library name matches the
            // user's wall clock. We don't need cross-machine sortability
            // here - the assay name is per-push-unique on a single
            // workstation, which a local timestamp already gives us.
            string assayName = $"Cadenza-{DateTime.Now:yyyyMMdd-HHmmss}";

            // Require the Skyline document to be saved so the BLIB
            // lives next to the .sky file at a stable on-disk path.
            // Writing to %TEMP% silently looked like it worked (the
            // RunCommand --add-library-* call reports success) but
            // Skyline's PeptideSettings.Libraries.LibrarySpecs change
            // did not stick - the Libraries list appeared empty in
            // the Peptide Settings dialog. Putting the BLIB next to
            // the .sky file at a permanent path avoids that and is
            // also where users expect their assay libraries to live.
            string? docPath = await TryGetSavedDocumentPathAsync();
            if (string.IsNullOrWhiteSpace(docPath))
            {
                StatusMessage = "Save the Skyline document first (File > Save), "
                    + "then re-run the push. Cadenza writes the assay BLIB next "
                    + "to the .sky file so Skyline can persist the library "
                    + "registration; an unsaved document doesn't have a stable "
                    + "location to anchor it to.";
                return;
            }
            blibPath = Path.Combine(Path.GetDirectoryName(docPath!)!, assayName + ".blib");

                // The schedule the Skyline-side outputs are built from. The
            // instrument method always uses the driver-only schedule:
            // passengers ride in windows that already exist, so including them
            // changes what the document extracts, never what the instrument does.
            var assaySchedule = IncludePassengers
                ? ScheduleResult.WithPassengers()
                : ScheduleResult;

            // Step 1: BLIB write + library registration.
            StatusMessage = $"Writing assay BLIB to {Path.GetFileName(blibPath)}...";
            var writeResult = await Task.Run(() =>
                BlibAssayWriter.Write(_candidates, assaySchedule, blibPath!, assayName));

            // AddSettingsListItem cannot register a Spectral Libraries
            // entry. Skyline's JsonToolServer.DeserializeSettingsItem
            // looks up `Deserialize(XmlReader)` on the generic argument
            // of SettingsListBase<T> - for the Spectral Libraries list
            // T is the abstract LibrarySpec, which has no static
            // Deserialize method (only its concrete subclasses do). So
            // we route through the SkylineCmd path instead via
            // RunCommand: --add-library-name + --add-library-path,
            // which calls CommandLine.SetLibrary(name, path). That
            // method invokes LibrarySpec.CreateFromPath internally,
            // which dispatches on the file extension to construct the
            // right concrete subclass (BiblioSpecLiteSpec for .blib)
            // and appends it to PeptideSettings.Libraries.LibrarySpecs
            // - the "active" list - in one shot.
            string registrationDiag;
            bool registered;
            try
            {
                string addLibOutput = await Task.Run(() =>
                    _skylineSession.Execute(c => c.RunCommand(new[]
                    {
                        // --add-library-path before --add-library-name:
                        // SkylineCmd is supposed to parse all args
                        // before running handlers (so the order
                        // shouldn't matter), but in practice RunCommand
                        // was leaving the library unattached when name
                        // came first. Trying path-first to see if the
                        // argument-stash order matters somewhere
                        // upstream.
                        "--add-library-path=" + blibPath!,
                        "--add-library-name=" + assayName,
                    })));
                if (string.IsNullOrWhiteSpace(addLibOutput))
                {
                    registered = true;
                    registrationDiag = "registered + activated via --add-library-name/--add-library-path";
                }
                else
                {
                    // Skyline's CommandLine writes failure reasons to
                    // stdout (file doesn't exist, name conflict, etc.).
                    // Any non-empty output means the command at minimum
                    // had something to say; we capture it verbatim and
                    // surface it - but don't bail, because some flags
                    // that succeed still echo informational text.
                    registered = !addLibOutput.Contains("Error", StringComparison.OrdinalIgnoreCase)
                              && !addLibOutput.Contains("conflict", StringComparison.OrdinalIgnoreCase);
                    registrationDiag = registered
                        ? $"--add-library-name/--add-library-path: {addLibOutput.Trim()}"
                        : $"library registration failed: {addLibOutput.Trim()}";
                }
            }
            catch (Exception ex)
            {
                registered = false;
                registrationDiag = $"RunCommand --add-library-* threw: {ex.Message}";
            }

            if (!registered)
            {
                StatusMessage = $"BLIB written to {blibPath} "
                    + $"({writeResult.RefSpectraWritten:n0} peptides), "
                    + $"but registering it with Skyline failed: {registrationDiag}. "
                    + "Add it manually via Settings > Peptide Settings > Library, "
                    + "then re-run the push so the target tree populates.";
                return;
            }

            // Step 1.5: align the document's peptide + transition
            // filter with the assay. Per Nick Shulman, RunCommand via
            // JSON-RPC behaves as if the user typed the SkylineCmd
            // flags into Skyline's Immediate Window, so the same
            // --tran-* / --peptide-* flags apply to the live document.
            // The most common cause of an empty target tree is a
            // default product-ion-charges = {1} filter, which silently
            // rejects every +2 fragment DIA-NN picked.
            var scheduledCands = assaySchedule.ScheduledIndices
                .Select(i => _candidates[i])
                .ToList();
            StatusMessage = $"BLIB registered ({registrationDiag}): "
                + $"{writeResult.RefSpectraWritten:n0} peptides. "
                + "Aligning document peptide + transition filter...";
            var settingsResult = await SkylineSettingsConfigurator.ConfigureAsync(
                _skylineSession!, scheduledCands, Mode,
                firingPadMin: FiringPadSec / 60.0);
            var recommendation = settingsResult.Recommendation;
            string settingsTail = string.IsNullOrWhiteSpace(settingsResult.SkylineOutput)
                ? "OK"
                : settingsResult.SkylineOutput;
            StatusMessage = $"Filter applied ({recommendation.ToStatusLine()}; Skyline: {settingsTail}). "
                + "Pushing transition list...";

            // Step 2: transition-list import. Each row carries the real
            // FragmentIon.Charge so Skyline's product-ion match succeeds
            // even for the +2 y/b fragments DIA-NN selects on longer
            // tryptic peptides.
            var csv = PeptideTransitionListBuilder.Build(_candidates, assaySchedule);
            int rowCount = csv.AsSpan().Count('\n') - 1;
            // Heavy isotope-labeled precursors (e.g. PRTC standards) are
            // omitted from the transition list by the builder - Skyline's
            // headless importer can't apply the document's explicit heavy
            // label and one such row aborts the whole import. They still
            // reach the assay as heavy through the registered library.
            int heavySkipped = scheduledCands.Count(c => !c.IsLightLabel);
            string heavyNote = heavySkipped > 0
                ? $" ({heavySkipped:n0} heavy-labeled precursors omitted from the list; "
                    + "they enter via the registered library)"
                : string.Empty;
            tlistPath = Path.Combine(Path.GetTempPath(),
                $"cadenza-targets-{Guid.NewGuid():N}.csv");
            await File.WriteAllTextAsync(tlistPath, csv);

            var pushStart = DateTimeOffset.UtcNow;
            const double rowsPerSec = 10.0;
            double etaMin = rowCount / rowsPerSec / 60.0;
            string etaText = etaMin < 1.0
                ? $"~{rowCount / rowsPerSec:0} s"
                : $"~{etaMin:0.0} min";
            StatusMessage = $"BLIB registered. Importing transition list: "
                + $"{assaySchedule.ScheduledIndices.Length:n0} precursors "
                + $"({rowCount:n0} rows){heavyNote}. ETA {etaText}. "
                + "Live progress in Skyline's Immediate Window.";

            using var heartbeatCts = new CancellationTokenSource();
            var heartbeat = HeartbeatStatus(pushStart, rowCount, rowsPerSec, heartbeatCts.Token);
            string tlistOutput = string.Empty;
            try
            {
                tlistOutput = await Task.Run(() => _skylineSession.Execute(c => c.RunCommand(new[]
                {
                    "--import-transition-list=" + tlistPath,
                    // Degrade any per-row rejection to a warning instead
                    // of aborting the whole import, so one unexpected row
                    // can't discard the entire scheduled assay.
                    "--ignore-transition-errors",
                })));
            }
            finally
            {
                heartbeatCts.Cancel();
                try { await heartbeat; } catch (OperationCanceledException) { }
            }

            string head = string.IsNullOrWhiteSpace(tlistOutput)
                ? "(no output)"
                : tlistOutput.Length > 400 ? tlistOutput.Substring(0, 400) + "..." : tlistOutput;
            var elapsed = DateTimeOffset.UtcNow - pushStart;
            StatusMessage = $"Pushed '{assayName}': "
                + $"BLIB at {Path.GetFileName(blibPath)} ({writeResult.RefSpectraWritten:n0} peptides) + "
                + $"transition list ({rowCount:n0} rows in {elapsed.TotalMinutes:0.0} min). "
                + $"Filter set to {recommendation.ToStatusLine()}. "
                + $"Skyline: {head}";
        }
        catch (Exception ex)
        {
            StatusMessage = $"Skyline push failed ({ex.GetType().Name}): {ex.Message}";
        }
        finally
        {
            if (tlistPath is not null)
            {
                try { File.Delete(tlistPath); } catch { /* best effort */ }
            }
        }
    }

    /// <summary>
    /// Updates <see cref="StatusMessage"/> every 10 s while a long Skyline
    /// transition-list import is in flight, so the user sees Cadenza is
    /// still alive while <c>RunCommand</c> blocks on the underlying
    /// <c>MassListImporter</c>.
    /// </summary>
    private async Task HeartbeatStatus(DateTimeOffset start, int rowCount, double rowsPerSec, CancellationToken ct)
    {
        try
        {
            while (!ct.IsCancellationRequested)
            {
                await Task.Delay(TimeSpan.FromSeconds(10), ct);
                var elapsed = DateTimeOffset.UtcNow - start;
                double pctEst = Math.Min(99, elapsed.TotalSeconds * rowsPerSec / Math.Max(1, rowCount) * 100);
                StatusMessage =
                    $"Skyline transition-list import in flight: elapsed {elapsed.TotalMinutes:0.0} min, "
                    + $"estimated {pctEst:0}% of {rowCount:n0} rows. "
                    + "Live per-protein progress in Skyline's Immediate Window.";
            }
        }
        catch (OperationCanceledException) { }
    }

    /// <summary>
    /// Query Skyline for the active document's on-disk path. Returns
    /// null if the document hasn't been saved (Skyline returns an
    /// empty path for unsaved documents), if the path doesn't actually
    /// exist on disk, or if the RPC call fails.
    /// </summary>
    private async Task<string?> TryGetSavedDocumentPathAsync()
    {
        try
        {
            string? docPath = await Task.Run(() => _skylineSession!.Execute(c => c.GetDocumentPath()));
            if (string.IsNullOrWhiteSpace(docPath)) return null;
            if (!File.Exists(docPath)) return null;
            return docPath;
        }
        catch
        {
            return null;
        }
    }


    [RelayCommand]
    private async Task LoadFromSkylineAsync()
    {
        if (_skylineSession is null)
        {
            StatusMessage = "Skyline: not connected. Launch this tool from Skyline's Tools menu.";
            return;
        }
        IsBusy = true;
        SkylineLoadCancellable = true;
        var pollWasRunning = _skylinePoll;
        _skylinePoll?.Dispose();
        _skylinePoll = null;
        _skylineLoadCts?.Dispose();
        _skylineLoadCts = new CancellationTokenSource();
        var ct = _skylineLoadCts.Token;

        try
        {
            var progress = new Progress<string>(msg =>
            {
                StatusMessage = msg;
                SkylineStatus = msg;
            });
            string? rtOverride = string.IsNullOrWhiteSpace(SkylineRtColumnOverride) ? null : SkylineRtColumnOverride;
            var loadResult = await SkylineLibraryLoader.LoadAsync(
                _skylineSession!, rtColumnOverride: rtOverride,
                qValueCutoff: QValueCutoff,
                progress: progress, cancellationToken: ct);
            if (loadResult.Candidates.Count == 0)
            {
                string empty = "Skyline: no candidates built from the document. "
                    + $"Probed RT column '{loadResult.RtColumnUsed}', fetched {loadResult.RawRowsFetched:n0} rows.";
                StatusMessage = empty;
                SkylineStatus = empty;
                return;
            }

            _candidates = loadResult.Candidates;
            _carafePrecursors = null; _carafeDetectability = null; // Skyline pool uses BLIB / measured peak widths
            _combinedDiannRows = null; _combinedFrags = null;
            _knownAccessions = _candidates.Select(c => c.ProteinGroup).ToHashSet();
            _geneToAccession = new Dictionary<string, IReadOnlyList<string>>();

            PerRunCurves = await Task.Run(() => CoverageCurves.PerRunCurves(_candidates));
            FullLibraryHeatmap = await Task.Run(() => CoverageCurves.BuildHeatmap(_candidates));

            LibraryPrecursors = _candidates.Count;
            LibraryPeptides = _candidates.Select(c => c.StrippedSequence).Distinct().Count();
            LibraryProteinGroups = _knownAccessions.Count;
            LibraryCarafeMatched = _candidates.Count(c => c.Fragments.Length > 0);

            LibraryLoaded?.Invoke();
            string rtSource;
            if (loadResult.BlibPathsConsulted.Count == 0)
            {
                rtSource = "all firing windows synthesised";
            }
            else
            {
                int total = loadResult.PeptidesFromBlibBoundaries + loadResult.PeptidesFromSynthesizedBoundaries;
                rtSource = $"{loadResult.PeptidesFromBlibBoundaries:n0} of {total:n0} firing windows from BLIB "
                    + $"({loadResult.BlibPathsConsulted.Count} libr{(loadResult.BlibPathsConsulted.Count == 1 ? "y" : "ies")})";
            }
            string ok = $"Skyline document loaded: {LibraryPrecursors:n0} precursors / {LibraryPeptides:n0} peptides / {LibraryProteinGroups:n0} groups (RT col '{loadResult.RtColumnUsed}', {loadResult.RawRowsFetched:n0} rows; {rtSource}).";
            StatusMessage = ok;
            SkylineStatus = ok;
            await RescheduleAsync();
            return;
        }
        catch (OperationCanceledException)
        {
            string cancelled = "Skyline: document load cancelled.";
            StatusMessage = cancelled;
            SkylineStatus = cancelled;
            return;
        }
        catch (Exception ex)
        {
            string err = $"Skyline document load failed: {ex.GetType().Name}: {ex.Message}";
            StatusMessage = err;
            SkylineStatus = err;
            return;
        }
        finally
        {
            IsBusy = false;
            SkylineLoadCancellable = false;
            _skylineLoadCts?.Dispose();
            _skylineLoadCts = null;
            if (pollWasRunning is not null)
            {
                _skylinePoll = new System.Threading.Timer(_ => PollSkyline(),
                    null, TimeSpan.FromSeconds(2), TimeSpan.FromSeconds(2));
            }
        }

        // Legacy probe path below for reference - unreachable because of the return above.
#pragma warning disable CS0162
        var unused = _skylinePoll;
        var pollWasRunningLegacy = _skylinePoll;
        try
        {
            // Probe a sequence of built-in report names with a per-call
            // timeout. Reports on a 72k-precursor document can take many
            // seconds to build server-side, so the timeout exists to keep
            // the UI responsive, not to give up on real work.
            string[] candidates =
            {
                "Peptide Transition List",
                "Transition Results",
                "Peptide Ratio Results",
                "Molecule Transition List",
            };
            var perProbeTimeout = TimeSpan.FromSeconds(30);

            SkylineTool.ReportRowsResult? firstHit = null;
            string? firstHitName = null;
            foreach (var name in candidates)
            {
                string msg = $"Skyline: fetching '{name}' (timeout {perProbeTimeout.TotalSeconds:0}s)...";
                System.Windows.Application.Current?.Dispatcher.Invoke(() =>
                {
                    StatusMessage = msg;
                    SkylineStatus = msg;
                });
                try
                {
                    var fetch = Task.Run(() => _skylineSession.Execute(c => c.GetReportRows(
                        reportName: name,
                        offset: 0,
                        count: 50,
                        columns: null,
                        filter: null,
                        includeMaxLength: false,
                        culture: "invariant")));
                    var winner = await Task.WhenAny(fetch, Task.Delay(perProbeTimeout));
                    if (winner != fetch)
                    {
                        string slow = $"Skyline: '{name}' didn't respond in {perProbeTimeout.TotalSeconds:0}s. Trying next.";
                        System.Windows.Application.Current?.Dispatcher.Invoke(() =>
                        {
                            StatusMessage = slow;
                            SkylineStatus = slow;
                        });
                        continue;
                    }
                    var result = await fetch;
                    if (result?.TotalRows > 0)
                    {
                        firstHit = result;
                        firstHitName = name;
                        break;
                    }
                    string empty = $"Skyline: '{name}' returned 0 rows. Trying next.";
                    System.Windows.Application.Current?.Dispatcher.Invoke(() =>
                    {
                        StatusMessage = empty;
                        SkylineStatus = empty;
                    });
                }
                catch (Exception probeEx)
                {
                    string rej = $"Skyline: '{name}' rejected ({probeEx.GetType().Name}: {probeEx.Message}). Trying next.";
                    System.Windows.Application.Current?.Dispatcher.Invoke(() =>
                    {
                        StatusMessage = rej;
                        SkylineStatus = rej;
                    });
                }
            }

            if (firstHit is null)
            {
                string fail = "Skyline: no built-in report returned data. Open Tools > Reports in Skyline "
                    + "and tell me the exact report name + columns I should map.";
                System.Windows.Application.Current?.Dispatcher.Invoke(() =>
                {
                    StatusMessage = fail;
                    SkylineStatus = fail;
                });
                return;
            }

            int cols = firstHit.Columns?.Length ?? 0;
            var colNames = firstHit.Columns is null
                ? "(no columns)"
                : string.Join(", ", firstHit.Columns.Select(c => c.Name));
            string ok = $"Skyline: '{firstHitName}' -> {firstHit.TotalRows:n0} rows / {cols} cols. "
                + $"Columns: {colNames}. "
                + "Send me this column list and the candidate builder will be wired against it.";
            System.Windows.Application.Current?.Dispatcher.Invoke(() =>
            {
                StatusMessage = ok;
                SkylineStatus = ok;
            });
        }
        catch (Exception ex)
        {
            string err = $"Skyline report fetch failed: {ex.GetType().Name}: {ex.Message}";
            System.Windows.Application.Current?.Dispatcher.Invoke(() =>
            {
                StatusMessage = err;
                SkylineStatus = err;
            });
        }
        finally
        {
            IsBusy = false;
            // Resume the background poll loop.
            if (pollWasRunningLegacy is not null)
            {
                _skylinePoll = new System.Threading.Timer(_ => PollSkyline(),
                    null, TimeSpan.FromSeconds(2), TimeSpan.FromSeconds(2));
            }
        }
#pragma warning restore CS0162
    }

    [RelayCommand]
    private void BrowseTargetList()
    {
        var dlg = new OpenFileDialog
        {
            Filter = "FASTA (*.fasta;*.fa)|*.fasta;*.fa|Text list (*.txt;*.csv)|*.txt;*.csv|All files (*.*)|*.*",
            Title = "Open target protein list (FASTA or accession/gene list)",
        };
        if (dlg.ShowDialog() == true)
        {
            try { TargetProteinsText = File.ReadAllText(dlg.FileName); }
            catch (Exception ex) { StatusMessage = $"Could not read {dlg.FileName}: {ex.Message}"; }
        }
    }

    [RelayCommand]
    private async Task LoadAsync()
    {
        bool hasReport = !string.IsNullOrEmpty(ReportPath) && File.Exists(ReportPath);
        bool hasCarafe = !string.IsNullOrEmpty(CarafeTsvPath) && File.Exists(CarafeTsvPath);
        bool hasTargets = !string.IsNullOrWhiteSpace(TargetProteinsText);

        if (!hasReport && !(hasCarafe && hasTargets))
            {
            StatusMessage = "Need either a DIA-NN report.parquet, OR a Carafe library + a target protein list.";
                return;
            }

        IsBusy = true;
        try
        {
            if (hasReport) await LoadFromDiannAsync();
            else await LoadFromCarafeOnlyAsync();
        }
        finally
        {
            IsBusy = false;
        }
    }

    private async Task LoadFromDiannAsync()
    {
        StatusMessage = "Reading DIA-NN report...";
        var rows = await DiannParquetReader.LoadAsync(ReportPath!, QValueCutoff);
        StatusMessage = $"Loaded {rows.Count:n0} precursors. Running parsimony...";

        bool hasCarafe = !string.IsNullOrEmpty(CarafeTsvPath) && File.Exists(CarafeTsvPath);
        _geneToAccession = ProteinListParser.BuildGeneToAccession(rows);

        // Carafe fragment overlay for the DIA-NN (experimental) rows. Present
        // on both the plain and the backfill path - it only swaps fragments.
        var frags = new Dictionary<FragmentKey, FragmentIon[]>();
        if (hasCarafe)
        {
            StatusMessage = $"Streaming Carafe library {Path.GetFileName(CarafeTsvPath)}...";
            var allowed = new HashSet<string>(rows.Select(r => CarafeKey.FromDiann(r.ModifiedSequence)));
            frags = await Task.Run(() => CarafeTsvReader.ExtractTopFragments(
                CarafeTsvPath!, allowed, Candidate.FragmentLimit));
        }

        if (BackfillWithPredictions && hasCarafe)
        {
            await LoadCombinedBackfillAsync(rows, frags);
            return;
        }

        var pepMap = DiannParquetReader.BuildPeptideProteinMap(rows);
        var parsimony = ParsimonyEngine.Assign(pepMap);
        _knownAccessions = parsimony.Values.Select(v => v.Group).ToHashSet();

        _candidates = CandidateBuilder.Build(rows, parsimony, frags);
        // Clear every Carafe-pool cache so the sliders no-op on this pool.
        _carafePrecursors = null; _carafeDetectability = null;
        _combinedDiannRows = null; _combinedFrags = null;
        await FinalizeLoadAsync(
            fragSourceLabel: frags.Count > 0 ? "Carafe library" : "DIA-NN report fragments",
            sourceLabel: "DIA-NN report");
    }

    /// <summary>
    /// Combined DIA-NN + Carafe "backfill" pool: confirmed DIA-NN targets
    /// (experimental tier) plus Carafe-predicted targets (predicted tier) that
    /// backfill leftover schedule capacity. The Carafe library is scoped to the
    /// proteins DIA-NN saw plus the optional target-protein list.
    /// </summary>
    private async Task LoadCombinedBackfillAsync(
        List<DiannRow> rows, Dictionary<FragmentKey, FragmentIon[]> frags)
    {
        // Scope predictions to proteins DIA-NN observed, plus any target list.
        var allowedAccessions = new HashSet<string>(StringComparer.Ordinal);
        foreach (var r in rows)
            foreach (var id in r.ProteinIds.Split(';', StringSplitOptions.RemoveEmptyEntries))
                allowedAccessions.Add(ProteinListParser.ExtractAccession(id.Trim()));
        foreach (var t in ProteinListParser.ParseText(TargetProteinsText))
            allowedAccessions.Add(t);

        StatusMessage = $"Streaming Carafe predictions for {allowedAccessions.Count:n0} proteins...";
        var carafePrecursors = await Task.Run(() =>
            CarafeLibraryReader.LoadCandidates(CarafeTsvPath!, allowedAccessions));

        // Parsimony is run over the DIA-NN evidence ONLY, so the confirmed
        // grouping and unique/razor classification are identical to a DIA-NN-only
        // run - the experimental schedule stays locked regardless of the Carafe
        // library. Predicted (Carafe) peptides are then routed by accession into
        // whichever confirmed group their protein occupies (or grouped among
        // themselves when their protein was unseen by DIA-NN). Merging predicted
        // entries first and confirmed second keeps the confirmed entries
        // authoritative for any peptide seen by both.
        var (expParsimony, expAccToGroup) = ParsimonyEngine.AssignWithAccessionGroups(
            DiannParquetReader.BuildPeptideProteinMap(rows));
        var predictedParsimony = CandidateBuilder.MapPredictedToExperimentalGroups(
            CandidateBuilder.CarafePeptideProteinMap(carafePrecursors), expAccToGroup);
        var parsimony = new Dictionary<string, (string Group, string Type)>(predictedParsimony);
        foreach (var (peptide, groupType) in expParsimony) parsimony[peptide] = groupType;
        _knownAccessions = parsimony.Values.Select(v => v.Group).ToHashSet();

        // Optional Bromo within-protein ranking overlay for the predicted tier.
        PredictedRankingReader.Result? predicted = null;
        if (!string.IsNullOrWhiteSpace(PredictedRankingPath))
        {
            StatusMessage = $"Reading Bromo ranking {Path.GetFileName(PredictedRankingPath)}...";
            predicted = await Task.Run(() =>
                PredictedRankingReader.LoadFromPath(PredictedRankingPath!));
        }

        // Optional Carafe detectability filter for the predicted tier.
        PredictedRankingReader.Result? detectability = null;
        var detPaths = ParsePredictedRankingPaths(DetectabilityPath);
        if (detPaths.Count > 0)
        {
            StatusMessage = detPaths.Count == 1
                ? $"Reading detectability {Path.GetFileName(detPaths[0])}..."
                : $"Reading detectability from {detPaths.Count} files...";
            detectability = await Task.Run(() => PredictedRankingReader.LoadFromPaths(detPaths));
        }

        // Cache both tiers' inputs so the peak-half-width / detectability
        // sliders rebuild the whole two-tier pool live (no library re-stream).
        _combinedDiannRows = rows;
        _combinedFrags = frags;
        _carafePrecursors = carafePrecursors;
        _carafeParsimony = parsimony;
        _carafePredicted = predicted;
        _carafeDetectability = detectability;

        RebuildCarafeCandidates();

        var built = _candidates!;
        int confirmed = built.Count(c => c.Source == CandidateSource.Experimental);
        int predictedCount = built.Count - confirmed;
        string fragSourceLabel =
            $"{confirmed:n0} confirmed + {predictedCount:n0} predicted; peak half-width {PeakHalfWidthSec:0} s";
        if (predicted is not null || detectability is not null)
        {
            // Do NOT overwrite the user's peptide-ranking choice here. In backfill
            // mode the confirmed tier must keep ranking by the user's evidence
            // metric (PrecursorQuantity / QValue) so the confirmed schedule is
            // identical to a backfill-off run; the scheduler ranks the PREDICTED
            // tier by its detectability / Bromo score independently. (Carafe-only
            // mode, where every peptide is predicted, still auto-selects
            // PredictedRanking in LoadFromCarafeOnlyAsync.)
            string src = predicted is not null ? "Bromo ranking" : "Carafe detectability ranking";
            fragSourceLabel += $"; predicted tier ranked by {src}";
        }
        if (detectability is not null)
            fragSourceLabel += $"; detectability >= {DetectabilityCutoff:0.00}";

        await FinalizeLoadAsync(
            fragSourceLabel: fragSourceLabel,
            sourceLabel: "DIA-NN report + Carafe backfill");
    }

    private async Task LoadFromCarafeOnlyAsync()
    {
        StatusMessage = "Parsing target protein list...";
        var tokens = ProteinListParser.ParseText(TargetProteinsText);
        var allowed = new HashSet<string>(tokens);
        if (allowed.Count == 0)
        {
            StatusMessage = "Target protein list is empty - cannot load Carafe-only mode.";
            return;
        }

        StatusMessage = $"Streaming Carafe library for {allowed.Count:n0} target accessions...";
        var carafePrecursors = await Task.Run(() =>
            CarafeLibraryReader.LoadCandidates(CarafeTsvPath!, allowed));
        if (carafePrecursors.Count == 0)
        {
            StatusMessage = "No Carafe precursors matched the target protein list. Check the accession format.";
            return;
        }

        StatusMessage = $"Loaded {carafePrecursors.Count:n0} Carafe precursors. Running parsimony...";
        var pepMap = CandidateBuilder.CarafePeptideProteinMap(carafePrecursors);
        var parsimony = ParsimonyEngine.Assign(pepMap);

        _geneToAccession = new Dictionary<string, IReadOnlyList<string>>();
        _knownAccessions = parsimony.Values.Select(v => v.Group).ToHashSet();

        // Optional predicted peptide-ranking overlay (Bromo / Carafe detectability).
        // One or more files and/or folders, one per line in the input box.
        PredictedRankingReader.Result? predicted = null;
        if (!string.IsNullOrWhiteSpace(PredictedRankingPath))
        {
            StatusMessage = $"Reading Bromo ranking {Path.GetFileName(PredictedRankingPath)}...";
            predicted = await Task.Run(() =>
                PredictedRankingReader.LoadFromPath(PredictedRankingPath!));
        }

        // Optional Carafe detectability filter overlay (separate input).
        PredictedRankingReader.Result? detectability = null;
        var detPaths = ParsePredictedRankingPaths(DetectabilityPath);
        if (detPaths.Count > 0)
        {
            StatusMessage = detPaths.Count == 1
                ? $"Reading detectability {Path.GetFileName(detPaths[0])}..."
                : $"Reading detectability from {detPaths.Count} files...";
            detectability = await Task.Run(() => PredictedRankingReader.LoadFromPaths(detPaths));
        }

        // Cache the raw inputs so the peak-half-width and detectability-cutoff
        // sliders can rebuild candidates + reschedule live without re-streaming.
        _carafePrecursors = carafePrecursors;
        _carafeParsimony = parsimony;
        _carafePredicted = predicted;
        _carafeDetectability = detectability;
        _combinedDiannRows = null; _combinedFrags = null; // Carafe-only: single tier

        RebuildCarafeCandidates();

        string fragSourceLabel = $"Carafe predicted fragments; peak half-width {PeakHalfWidthSec:0} s";
        // Auto-select the predicted ranking whenever either overlay supplies a
        // within-protein priority: Bromo takes precedence, else detectability.
        if (predicted is not null || detectability is not null)
        {
            int matched = _candidates!.Count(c => c.PredictedPriority is not null);
            PeptideRanking = PeptidePriority.PredictedRanking;
            string src = predicted is not null ? "Bromo ranking" : "Carafe detectability ranking";
            fragSourceLabel = $"{src} matched {matched:n0}/{_candidates!.Count:n0}";
        }
        if (detectability is not null)
            fragSourceLabel += $"; detectability >= {DetectabilityCutoff:0.00} ({_candidates!.Count:n0} precursors kept)";

        await FinalizeLoadAsync(
            fragSourceLabel: fragSourceLabel,
            sourceLabel: "Carafe library + target list");
    }

    /// <summary>
    /// Rebuild the Carafe candidate pool from the cached precursors using the
    /// current peak-half-width, predicted-ranking overlay, and detectability
    /// filter. Cheap (no library re-stream); used at load and by the sliders.
    /// </summary>
    private void RebuildCarafeCandidates()
    {
        if (_carafePrecursors is null || _carafeParsimony is null) return;
        if (_combinedDiannRows is not null && _combinedFrags is not null)
        {
            // Combined backfill pool: experimental DIA-NN tier + predicted
            // Carafe tier, deduped on (peptide, charge).
            _candidates = CandidateBuilder.BuildCombined(
                _combinedDiannRows, _combinedFrags, _carafePrecursors, _carafeParsimony,
                PeakHalfWidthMin, _carafePredicted, _carafeDetectability, DetectabilityCutoff);
            return;
        }
        _candidates = CandidateBuilder.BuildFromCarafe(
            _carafePrecursors, _carafeParsimony, PeakHalfWidthMin,
            _carafePredicted, _carafeDetectability, DetectabilityCutoff);
    }

    /// <summary>Rebuild the Carafe pool and refresh stats + schedule (slider path).</summary>
    private async Task RefreshCarafeAsync()
    {
        if (_carafePrecursors is null) return;
        RebuildCarafeCandidates();
        LibraryPrecursors = _candidates!.Count;
        LibraryPeptides = _candidates.Select(c => c.StrippedSequence).Distinct().Count();
        LibraryProteinGroups = _candidates.Select(c => c.ProteinGroup).Distinct().Count();
        await RescheduleAsync();
    }

    private async Task FinalizeLoadAsync(string fragSourceLabel, string sourceLabel)
    {
        if (_candidates is null) return;

        PerRunCurves = await Task.Run(() => CoverageCurves.PerRunCurves(_candidates));
        FullLibraryHeatmap = await Task.Run(() => CoverageCurves.BuildHeatmap(_candidates));

        LibraryPrecursors = _candidates.Count;
        LibraryPeptides = _candidates.Select(c => c.StrippedSequence).Distinct().Count();
        LibraryProteinGroups = _knownAccessions?.Count ?? 0;
        LibraryCarafeMatched = _candidates.Count(c => c.Fragments.Length > 0);

        LibraryLoaded?.Invoke();
        StatusMessage = $"Source: {sourceLabel}. Loaded {LibraryPrecursors:n0} precursors / {LibraryPeptides:n0} peptides / {LibraryProteinGroups:n0} groups. Fragments: {fragSourceLabel}.";
        await RescheduleAsync();
    }

    [RelayCommand]
    public async Task RescheduleAsync()
    {
        if (_candidates is null) return;
        _rescheduleCts?.Cancel();
        var cts = _rescheduleCts = new CancellationTokenSource();
        var p = CurrentParameters();
        try
        {
            var result = await Task.Run(() => Scheduler.Run(_candidates!, p, cts.Token), cts.Token);
            if (cts.IsCancellationRequested) return;
            ApplyResult(result, p);
        }
        catch (System.OperationCanceledException) { }
    }

    private void ApplyResult(ScheduleResult result, SchedulingParameters p)
    {
        ScheduleResult = result;
        ScheduledPrecursors = result.ScheduledIndices.Length;
        SlotsTotal = result.Slots.Length;
        SlotsShared = result.Slots.Count(s => s.MemberIndices.Count > 1);
        PeakLoad = result.SlotCountCurve.Length == 0 ? 0 : result.SlotCountCurve.Max();
        ProteinsCovered = result.ProteinGroupsCovered;

        // Per-group rollups for stats + coverage plot.
        var scheduledByGroup = new Dictionary<string, List<int>>();
        for (int i = 0; i < result.ScheduledIndices.Length; i++)
        {
            int candIdx = result.ScheduledIndices[i];
            var g = _candidates![candIdx].ProteinGroup;
            if (!scheduledByGroup.TryGetValue(g, out var list))
            {
                list = new List<int>();
                scheduledByGroup[g] = list;
            }
            list.Add(candIdx);
        }

        // Distinct peptide-sequence count across scheduled precursors.
        var seenPep = new HashSet<string>();
        foreach (var idx in result.ScheduledIndices)
            seenPep.Add(_candidates![idx].StrippedSequence);
        ScheduledPeptides = seenPep.Count;

        // Median peptides per scheduled group (excluding 0-coverage).
        var perGroup = scheduledByGroup.Select(kv => kv.Value.Select(i => _candidates![i].StrippedSequence).Distinct().Count()).ToList();
        perGroup.Sort();
        MedianPeptidesPerProtein = perGroup.Count == 0 ? 0
            : perGroup.Count % 2 == 1 ? perGroup[perGroup.Count / 2]
            : (perGroup[perGroup.Count / 2 - 1] + perGroup[perGroup.Count / 2]) / 2.0;

        // Coverage plot data: one row per protein group with summed
        // peptide intensity across the entire library + peptides actually
        // scheduled. Sort descending by intensity for easy display.
        var groupIntensity = new Dictionary<string, double>();
        foreach (var c in _candidates!)
        {
            if (!groupIntensity.TryGetValue(c.ProteinGroup, out var sum)) sum = 0.0;
            groupIntensity[c.ProteinGroup] = sum + c.PrecursorQuantity;
        }
        ProteinCoverage.Clear();
        ProteinCoverage.AddRange(
            groupIntensity
                .OrderByDescending(kv => kv.Value)
                .Select(kv => new ProteinCoverageRow(
                    kv.Key,
                    kv.Value,
                    scheduledByGroup.TryGetValue(kv.Key, out var idxs)
                        ? idxs.Select(i => _candidates![i].StrippedSequence).Distinct().Count()
                        : 0)));

        // Build the scheduled-only heatmap so the user can toggle between
        // "full library" and "what's actually in the assay".
        if (_candidates is { Count: > 0 } && result.ScheduledIndices.Length > 0)
        {
            var subset = new List<Candidate>(result.ScheduledIndices.Length);
            foreach (var idx in result.ScheduledIndices) subset.Add(_candidates[idx]);
            ScheduledHeatmap = CoverageCurves.BuildHeatmap(subset);
        }
        else
        {
            ScheduledHeatmap = null;
        }

        StatusMessage = $"Scheduled {ScheduledPrecursors:n0} precursors / {ScheduledPeptides:n0} peptides covering {ProteinsCovered:n0}/{_knownAccessions?.Count:n0} groups in {SlotsTotal:n0} slots (peak load {PeakLoad}/{p.CycleBudget}).";
        ScheduleUpdated?.Invoke();
    }

    partial void OnHeatmapViewChanged(HeatmapView value)
    {
        OnPropertyChanged(nameof(Heatmap));
        HeatmapViewChanged?.Invoke();
    }

    /// <summary>Fires when the heatmap view toggle changes, so the View can re-render.</summary>
    public event Action? HeatmapViewChanged;

    private SchedulingParameters CurrentParameters()
    {
        var targets = new HashSet<string>();
        if (!string.IsNullOrWhiteSpace(TargetProteinsText) && _knownAccessions is not null && _geneToAccession is not null)
        {
            var tokens = ProteinListParser.ParseText(TargetProteinsText);
            var resolved = ProteinListParser.ResolveGenes(tokens, _geneToAccession, _knownAccessions);
            foreach (var t in resolved) targets.Add(t);
        }
        TargetProteinsParsedCount = targets.Count;
        return new SchedulingParameters
        {
            Mode = Mode,
            Objective = Objective,
            IsolationWindowTh = IsolationWindowTh,
            FragmentTolDa = FragmentTolDa,
            CycleBudget = CycleBudget,
            FiringPadSec = FiringPadSec,
            QValueCutoff = QValueCutoff,
            PrmIsolationWidthTh = PrmIsolationWidthTh,
            ChargeHandling = ChargeHandling,
            NormalizedCollisionEnergy = NormalizedCollisionEnergy,
            TargetProteins = targets,
            TargetMode = TargetMode,
            MinPeptidesPerProtein = MinPeptidesPerProtein,
            MaxPeptidesPerProtein = MaxPeptidesPerProtein,
            PeptideRankCap = PeptideRankCap,
            MinTargetIonShare = MinTargetIonShare,
            // Always on: every top-4 fragment of a target must carry under 1%
            // co-isolated signal from identified precursors, background included.
            BackgroundFragmentClash = true,
            ProteinRanking = ProteinRanking,
            PeptideRanking = PeptideRanking,
        };
    }

    // Parameter-change hooks: each one debounces by simply firing a fresh
    // reschedule. RescheduleAsync cancels prior runs, so the slider can
    // move freely without queuing stale work.
    partial void OnModeChanged(AcquisitionMode value) => _ = RescheduleAsync();
    partial void OnObjectiveChanged(CoverageObjective value) => _ = RescheduleAsync();
    partial void OnIsolationWindowThChanged(double value) => _ = RescheduleAsync();
    partial void OnFragmentTolDaChanged(double value) => _ = RescheduleAsync();
    partial void OnCycleBudgetChanged(int value) => _ = RescheduleAsync();
    partial void OnFiringPadSecChanged(double value) => _ = RescheduleAsync();

    // Peak half-width (baked into RtStart/RtStop) and the detectability cutoff
    // (a candidate filter) both change the Carafe pool at build time, so they
    // rebuild candidates from the cached precursors rather than just
    // rescheduling. No-op for DIA-NN / Skyline pools (cache is null).
    partial void OnPeakHalfWidthSecChanged(double value)
    {
        OnPropertyChanged(nameof(PeakHalfWidthMin));
        _ = RefreshCarafeAsync();
    }

    partial void OnDetectabilityCutoffChanged(double value) => _ = RefreshCarafeAsync();

    partial void OnPrmIsolationWidthThChanged(double value) => _ = RescheduleAsync();
    partial void OnMinTargetIonShareChanged(double value) => _ = RescheduleAsync();
    partial void OnChargeHandlingChanged(ChargeHandling value) => _ = RescheduleAsync();
    // NormalizedCollisionEnergy is output-only - no reschedule.
    partial void OnMinPeptidesPerProteinChanged(int value) => _ = RescheduleAsync();
    partial void OnMaxPeptidesPerProteinChanged(int value) => _ = RescheduleAsync();
    partial void OnPeptideRankCapChanged(int value) => _ = RescheduleAsync();
    partial void OnProteinRankingChanged(ProteinPriority value) => _ = RescheduleAsync();
    partial void OnPeptideRankingChanged(PeptidePriority value) => _ = RescheduleAsync();
    partial void OnTargetModeChanged(TargetListMode value) => _ = RescheduleAsync();
    partial void OnTargetProteinsTextChanged(string value) => _ = RescheduleAsync();
}

/// <summary>One row in the protein-coverage plot's data series.</summary>
public sealed record ProteinCoverageRow(string ProteinGroup, double SummedIntensity, int PeptidesScheduled);

/// <summary>Which precursor set drives the m/z x RT heatmap.</summary>
public enum HeatmapView
{
    /// <summary>Every precursor in the loaded DIA-NN report.</summary>
    FullLibrary,
    /// <summary>Only precursors selected by the current schedule.</summary>
    ScheduledAssay,
}
