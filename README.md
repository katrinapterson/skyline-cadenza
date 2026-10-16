<p align="center">
  <img src="docs/screenshots/logo.jpg" alt="Skyline Cadenza" width="540"/>
</p>

# Skyline Cadenza

A Skyline external tool for designing **load-balanced PRM / MTM targeted assays** from gas-phase-fractionated DIA-NN libraries or Carafe AI-predicted spectral libraries.

Given a DIA-NN `report.parquet` (single-run, GPF, or quantitative experiment), a Carafe library plus a target protein list, both together, or the active Skyline document, Cadenza builds a candidate pool, applies protein parsimony, and selects a peptide target list that **maximizes the number of parsimonious proteins covered** under a configurable concurrent-MS/MS-slot budget. Multiplexed Multiple Target Monitoring (MTM) slot-sharing follows the scheme in US 11,688,595 B2 (Remes & MacCoss, 2023): co-eluting precursors that fit in a 2-3 Th isolation window share one acquisition slot when their top-4 fragments do not overlap. On top of that, Cadenza judges every window by everything the quadrupole would transmit, targeted or not: a target is not placed where co-isolated ions would take too much of the fill or interfere with its fragments, and windows firing at the same time never overlap in m/z.

The tool exposes a live-update WPF UI with every parameter bound to a control, recomputes the schedule in the background on each change, and writes results back to the active Skyline document **or** to a Thermo Method Editor-importable inclusion list.

Analogous to Thermo's PRM-Conductor but focused on the MacCoss lab's MTM work; descended from the prototype notebook at [`maccoss/targeted-modeling`](https://github.com/maccoss/targeted-modeling) (`gpf_coverage.ipynb`).

## What it looks like

Cadenza loads on the left, with every scheduling knob bound to a live control. Plots are tabbed on the right; the **summary strip** along the top of the plot area tracks library counts vs. scheduled counts in real time:

<p align="center">
  <img src="docs/screenshots/cadenza-ui.png" alt="Skyline Cadenza main UI" width="900"/>
</p>

## What it produces

Each protein group is plotted by `log10(summed precursor intensity)` against rank and colored by the number of peptides the scheduler kept for it. High-intensity groups fill the scheduler easily (deep green = 5+ peptides each); low-intensity groups drop out (red / grey tail) when the cycle budget can't fit more. The title reports how many groups are covered:

<p align="center">
  <img src="docs/screenshots/cadenza-protein-coverage.png" alt="Protein coverage by peptides scheduled" width="900"/>
</p>

The **m/z × RT heatmap** in scheduled mode draws every slot as its emitted isolation window: width = member m/z span + the PRM / solo-slot width, height = padded firing window. Color encodes multiplex depth: blue = 1 target / slot (solo), then green, orange, red, and purple for 5+. This view is the fastest way to see how the scheduler multiplexes across the gradient:

<p align="center">
  <img src="docs/screenshots/cadenza-mtm-isolation.png" alt="MTM isolation windows on m/z x RT heatmap" width="900"/>
</p>

## Installation

1. Download `SkylineCadenza-v<version>.zip` from the [GitHub Releases](https://github.com/maccoss/skyline-cadenza/releases) page, or build it from source (see *Development* below; output is `publish/SkylineCadenza.zip`).
2. In Skyline-daily: `Tools > External Tools > Add > File...` and pick the zip.
3. The tool appears as `Skyline Cadenza...` in the Tools menu.

Requires **Skyline-daily 26.1.1.083 or later** and the **.NET 8 Desktop Runtime** on Windows.

## Quick start: DIA-NN flow

1. Launch via `Tools > Skyline Cadenza...`. The Skyline-status line should immediately show your current document and counts.
2. **Inputs panel** (left side):
   - **DIA-NN report.parquet**: the report from your library search.
   - **Carafe library TSV** *(optional)*: confirmed peptides always use the fragments DIA-NN measured (`Fr.0.Id`..`Fr.11.Id`). Load a Carafe library and tick **Backfill unused schedule time with Carafe predictions** to fill leftover cycle budget with predicted peptides (see *DIA-NN + Carafe backfill* below).
3. Click **Load + schedule**. The progress bar runs and the status line tells you which stage is in flight (read → parsimony → fragments → candidates → scheduler).
4. Tune the parameters live in the left panel. The scheduler reruns automatically after each change and the plot tabs redraw.
5. When you're happy with the schedule:
   - **Push targets to Skyline** (save the Skyline document first) writes an assay spectral library next to the document, registers it, aligns the document's peptide and transition settings with the assay, and imports the targets into the protein/peptide tree. See *Outputs* below.
   - **Export Thermo CSV...** writes a Thermo Method Editor-importable scheduled inclusion list. Open your template `.meth` in Method Editor and choose `Import Scheduled List`.

To schedule from peptides already in Skyline instead, click **Load from Skyline document**. Cadenza reads the document's precursors and takes each peptide's firing window from the measured peak boundaries in the document's BiblioSpec libraries, falling back to RT ± 18 s for peptides no library covers. **RT column override** picks the retention-time column explicitly when the automatic probe chooses the wrong one.

## Carafe-only mode (no DIA-NN report)

Leave the report path blank and instead provide:

- A **Carafe library TSV** (`carafe_spectral_library.tsv` or similar).
- A **target protein list**: paste a FASTA, accessions, or gene symbols into the text box, or click *Open FASTA / list...*.

Cadenza streams the Carafe TSV, filters to the target proteins by `ProteinID`, computes parsimony from the Carafe peptide-protein associations, and synthesizes peak boundaries from `Tr_recalibrated` ± the **Predicted peak half-width** (default 6 s).

Optional prediction inputs:

- **Bromo peptide ranking**: a `peptide_rankings.tsv` that sets the within-protein peptide priority.
- **Carafe detectability predictions**: one or more `*_detectability_pred.tsv` files or a folder of shards. Precursors below the **Detectability cutoff** (default 0.5) are dropped; peptides with no prediction are kept. With no Bromo file, detectability also ranks peptides within each protein.

Predictions carry no cross-protein abundance, so they only rank peptides within a protein, never proteins.

## DIA-NN + Carafe backfill

Provide both a DIA-NN report and a Carafe library and tick **Backfill unused schedule time with Carafe predictions**. Confirmed DIA-NN targets are scheduled first; Carafe-predicted peptides then fill whatever cycle budget is left, deduplicated against the confirmed targets by peptide and charge and limited to the proteins DIA-NN observed plus any target list. The confirmed schedule is identical whether or not a Carafe library is loaded: predictions only ever add targets.

## Parameter reference

All knobs are in the left panel and bound to live controls. The full reference, with the reasoning behind each default, is in [`docs/settings.md`](docs/settings.md).

### Target protein list
- Paste a FASTA, newline-separated accessions, or gene symbols, or *Open FASTA / list...* a file.
- **Target list first, fill with rest** (default): the scheduler covers your target proteins first, then fills remaining budget with off-list proteins.
- **Exclusive (target list only)**: ignores off-list proteins entirely.

### Forced peptides
- One peptide per line (`PEPTIDEK`, `PEPTIDEK+2`, or a charge column), or *Open peptide list...*. Scheduled before everything else, for carrying the same peptides across schedules. Forced peptides obey every quality and instrument rule but ignore ranking, the rank cap, the per-protein bounds, and an exclusive target list. An entry with no charge tries its best charge state first.

### Acquisition mode
- **MTM (multiplex)** (default): co-eluting precursors share isolation windows when the rules below allow.
- **PRM (1 precursor / slot)**: solo windows only.

### Coverage objective
- **Balanced**: the published Stellar-webinar algorithm. Takes each protein's best peptide that fits the budget, then adds more per protein round-robin up to *Maximum peptides per protein*.
- **Maximize Protein Coverage** (default): look-ahead cover that prefers peptides able to join an existing window for free, then the same round-robin load-up.
- **Maximize Peptide Coverage**: look-ahead cover with no per-protein limit on the load-up.

See [`docs/coverage-objectives.md`](docs/coverage-objectives.md).

### MTM constraints
- **Isolation window cap (Th)**: default 3.0; the widest window an MTM slot may grow to.
- **Fragment tolerance (Da)**: default 0.5; precursors with any top-4 fragment within this tolerance won't multiplex.
- **PRM / solo-slot window width (Th)**: default 0.7; the quadrupole width around each precursor. Every member's window must fit inside the slot, so the emitted width is the member m/z span plus this.
- **Minimum target share of window current**: default 1.5%, 0 = off. Refuses a placement that would leave any target below this share of the ion current its window admits while it elutes, counting every identified precursor in the band, targeted or not.
- **Max targets per slot**: default 0 (no cap). Caps how many targeted precursors one window may multiplex. Untargeted precursors and passengers don't count.
- **Charge handling per slot**:
  - *Match charge per slot* (default): an MTM slot can only contain precursors of identical charge.
  - *Allow mixed (use slot majority charge)*: lets different charges co-multiplex; the Thermo CSV reports the slot's majority charge.

Two rules are always on and have no control: a target's top-4 fragments must also be clear of co-isolated background precursors, and two windows firing at the same time never overlap in m/z (see [`docs/scheduler-window-scoping.md`](docs/scheduler-window-scoping.md)).

### Cycle budget
- **Concurrent slots**: default 100. Maximum simultaneous acquisition slots at any retention time.
- **RT buffer (each side)**: default 15 s. Padding added on each side of every peak's boundaries. A slot's firing window is the padded union of its members.

### Per-protein peptide bounds
- **Minimum peptides per protein**: default 1. Protein groups that end up below the minimum are dropped from the schedule and their slot capacity is freed. Groups holding a forced peptide are exempt.
- **Maximum peptides per protein**: default 5.

### Ranking
- **Protein ranking**: `SummedIntensity` (default), `ProteinQValue`, `ProvidedListOrder`.
- **Peptide ranking**: `PrecursorQuantity` (default), `QValue`, or the predicted ranking (selected automatically in Carafe-only mode when a Bromo or detectability file is loaded).
- **Rank cap: top N per protein**: default 0 (no cap). Only each protein's top N precursors by the peptide ranking are eligible at all.

### Library filter
- **Q value cutoff**: default 0.01. Applied at ingest before parsimony.

### Output
- **Normalized CE**: default 28.0. Only written into the Thermo CSV; doesn't affect scheduling.
- **Include passengers in Skyline push**: default off. See *Passengers* below.

## Plot tabs

| Tab | Library view shows | Scheduled view shows |
| --- | --- | --- |
| **Cycle load** | Concurrent acquisition slots vs. retention time, with budget line. | (single view) |
| **Protein coverage** | Per-protein bar of `log10(summed precursor intensity)`, sorted high → low, colored by *peptides scheduled* per protein (grey = 0, red = 1, ..., dark green = 5+). | (single view) |
| **Per-run coverage** | Coverage curves per DIA-NN run. If runs match the `Chrlib400-500` GPF pattern, one curve per GPF + dashed black sum. Otherwise one curve labeled with the run name. | (single view) |
| **m/z × RT heatmap** | Binned precursor density across the loaded library (**Full library density**). | **Scheduled isolation windows**: each slot drawn as its emitted window, width = member m/z span + PRM width, height = padded firing window. Color by target count (PRM = solid blue; MTM = blue/green/orange/red/purple for 1/2/3/4/5+). |
| **Slot occupancy** | (n/a) | Histogram of slots by number of targets, with max / median / mode in the title. |

## Outputs

### Push targets to Skyline

**Push targets to Skyline** requires a saved Skyline document and runs four steps over Skyline's JSON-RPC connection:

1. Writes a self-contained BiblioSpec `.blib` next to the document: one reference spectrum per scheduled precursor (top 6 fragments) with its firing window.
2. Registers that library in the document.
3. Aligns the document's peptide, transition, and full-scan settings with the assay (precursor and product charges, ion types, library ion picking, peptide length bounds, scheduling-window RT filter, MS1 isotopes, and DIA or PRM acquisition). Any flag the installed Skyline doesn't recognise is reported in the status line without aborting the push.
4. Imports a peptide-style transition list into the live document. Each precursor gets M+0, M+1, and M+2 precursor rows followed by up to six fragment rows with their real product charges. Columns:

```
Protein Name, Peptide Modified Sequence, Precursor Charge, Precursor m/z, Product m/z, Product Charge, Explicit Retention Time, Explicit Retention Time Window, Note
```

The `Note` column carries `slot=<id>;type=<unique|razor>` (plus `;passenger` for passengers) so MTM grouping is visible inside Skyline. DIA-NN's `C(UniMod:4)` and Carafe's `_C[UniMod:4]_` modified-sequence syntaxes are both normalised to Skyline's `C[UniMod:4]`. Heavy isotope-labeled standards are left out of the transition list; they reach the document through the registered library.

Two settings have no command-line flag and need a one-time manual step in **Peptide Settings > Prediction**: set the alignment target to library RTs, and check *Impute boundaries for missing peaks*. Both persist across Skyline sessions.

### Passengers

Passengers are precursors that already fall inside a scheduled window's isolation band and elute entirely within its firing time, with at least 5% of the band's co-isolated current. The instrument fragments them either way, so they cost nothing and never change the acquisition method. Tick **Include passengers in Skyline push** to add them to the Skyline document as well.

### Export Thermo CSV

The **Export Thermo CSV...** button writes a scheduled inclusion list for Thermo Method Editor. Columns:

```
Compound, Formula, Adduct, m/z, z, t start (min), t stop (min), Isolation Window (m/z), HCD Collision Energy
```

For **PRM** mode: one row per scheduled precursor. For **MTM** mode: **one row per slot**. Multiplexed slots join member peptide names in the `Compound` column with ` | ` separators and report the slot's center m/z, the majority charge in `z`, the padded firing window in `t start` / `t stop`, and the member m/z span plus the PRM width in `Isolation Window (m/z)`. Passengers are never written; they ride in the targets' windows.

Import the CSV into a template `.meth` via Thermo Method Editor > `Import Scheduled List`.

### Driving Cadenza from another app

A pipeline can get the same CSV without Skyline, either headless or by opening the GUI preloaded:

```
SkylineCadenzaCli export-thermo REPORT.parquet PARAMS.json|- OUTPUT.csv
SkylineCadenza.exe --report <report.parquet> [--export-dir <folder>]
```

`export-thermo` runs the app's DIA-NN Load and Thermo export with the app's default settings (`-`), or with the settings named in PARAMS.json. It prints a one-line summary to stderr. `--report` opens the GUI standalone, with no Skyline connection and Push to Skyline disabled. It loads and schedules the report on startup. `--export-dir` makes **Export Thermo CSV...** open in that folder. See [docs/automation.md](docs/automation.md) for parameters, defaults, and exit codes.

## Architecture

```
src/
  SkylineCadenza.App/      WPF UI (net8.0-windows). The .exe Skyline launches.
  SkylineCadenza.Core/     Pure algorithms + I/O + RPC client (net8.0-windows).
  SkylineCadenza.Cli/      Headless JSON-in / JSON-out scheduler for scripts
                           and the notebooks.
  SkylineCadenza.Tests/    xUnit.
external/SkylineTool/      Vendored copies of pwiz's JSON-RPC client sources
                           (link-compiled into Core; Apache 2.0).
docs/                      Settings, ingest, Skyline integration, and
                           scheduler design notes; screenshots.
notebooks/                 Algorithm-comparison notebook (uses the CLI).
release-notes/             Per-version release notes.
build/                     MSBuild packaging target -> SkylineCadenza.zip.
tools/                     Dev helpers (install-dev.ps1).
```

Highlights:
- **Parsimony**: streamlined port of `skyline_prism.parsimony` (inverted-index subsumption + single-pass razor assignment), equivalent to the reference on the validated fixtures.
- **Scheduler**: greedy two-pass (cover → load-up), strict intersection-based co-elution (no chain extension across non-overlapping peptides), explicit MTM slot model with co-isolation checks against the whole candidate pool. Deterministic: the same inputs and settings always give the same schedule.
- **Ingest**: `Parquet.Net` for DIA-NN parquet; hand-rolled streaming TSV reader for large Carafe libraries; direct SQLite reads of Skyline's BiblioSpec libraries.
- **CLI**: `SkylineCadenzaCli schedule INPUT.json OUTPUT.json` runs the real scheduler on a JSON request; `schedule-diann REPORT.parquet PARAMS.json OUTPUT.json` runs the app's DIA-NN ingest first; `export-thermo REPORT.parquet PARAMS.json|- OUTPUT.csv` writes the app's Thermo CSV (see [docs/automation.md](docs/automation.md)).
- **Skyline RPC**: connect-per-call against the JSON-RPC pipe (`SkylineMcpJson-` prefix) Skyline's `JsonToolServer` exposes, the same pattern the upstream `SkylineMcpServer` uses.
- **Plots**: ScottPlot 5 WPF.

## Development

Requires .NET 8 SDK on Windows or WSL (cross-compile is enabled in `Directory.Build.props`).

```bash
git clone https://github.com/maccoss/skyline-cadenza.git
cd skyline-cadenza
dotnet restore
dotnet build SkylineCadenza.sln
dotnet test
dotnet msbuild build/package.proj    # produces publish/SkylineCadenza.zip
```

For local smoke testing without going through the zip+install cycle, run `tools/install-dev.ps1` from Windows. It publishes the App and copies the output into `%LOCALAPPDATA%\Apps\SkylineDaily\Tools\SkylineCadenza\`. Restart Skyline-daily afterwards.

Releases are cut by pushing a `v<version>` tag; see [`release-notes/README.md`](release-notes/README.md).

## Status

| Area | Status |
| --- | --- |
| DIA-NN report ingest | ✅ stable |
| Carafe TSV ingest + detectability / Bromo ranking | ✅ stable |
| Carafe-only mode (target list + Carafe) | ✅ working |
| DIA-NN + Carafe backfill | ✅ working |
| Load from running Skyline document | ✅ working; firing windows from BLIB peak boundaries, falls back to RT ± 18 s |
| Parsimony engine | ✅ stable |
| Scheduler (PRM + MTM, co-isolation checks, forced peptides, passengers) | ✅ stable |
| Live UI + plots | ✅ working |
| Thermo CSV export | ✅ working |
| Push targets to Skyline (library + settings + target tree) | ✅ working; two Prediction settings still need a one-time manual step |
| Per-GPF heatmap panels | ⏳ single combined heatmap today; per-GPF stacked panels still pending |

## References

- Heil L. et al., *Closing the gap between targeted and untargeted measurements using intelligent data acquisition on Stellar MS*, ASMS 2025.
- Remes P.R. & MacCoss M.J., US 11,688,595 B2, *Mass spectrometry methods*, 2023.
- Stellar MS Webinar 2024, *Balancing the Load: maximizing instrument time and performance characteristic to increase targeted data quality*.
- ProteoWizard pwiz: <https://github.com/ProteoWizard/pwiz>.
- Prototype notebook: <https://github.com/maccoss/targeted-modeling> (`gpf_coverage.ipynb`).

## License

Apache 2.0. See `LICENSE`.
