# Driving Cadenza from another app

Two entry points let a pipeline (for example AutoTargeting, which runs
GPF-DIA discovery, DIA-NN, assay design, and the targeted run on a Stellar)
use Cadenza's assay design without Skyline:

| Situation | Entry point | Result |
|---|---|---|
| Unattended run | `SkylineCadenzaCli export-thermo` | Thermo CSV written directly |
| A person wants to tune the assay | `SkylineCadenza.exe --report ... --export-dir ...` | GUI opens with the report loaded and scheduled; the person exports into the watched folder |

Both produce the same Thermo Method Editor scheduled inclusion list as the
app's **Export Thermo CSV...** button (see the README's "Export Thermo CSV"
section for the columns). With the same report and settings the two are
byte-identical, because they run the same Core code.

## Headless: `export-thermo`

```
SkylineCadenzaCli export-thermo REPORT.parquet PARAMS.json OUTPUT.csv
SkylineCadenzaCli export-thermo REPORT.parquet - OUTPUT.csv
```

- `REPORT.parquet`: a DIA-NN `report.parquet`.
- `PARAMS.json`: scheduling settings, or `-` to use the app's defaults for
  every setting. Any field left out of the JSON also takes the app's default.
- `OUTPUT.csv`: the Thermo CSV. Missing parent folders are created.

The pipeline is the app's **Empirical (DIA-NN report)** Load followed by
**Export Thermo CSV...**. It reads the report at the q-value cutoff, runs
protein parsimony, builds candidates from the report's own fragments,
resolves the optional target protein list, schedules, and writes the CSV
(`DiannThermoExport.RunAsync` in Core). No Carafe library or Skyline
document is involved.

On success it prints one line to stderr and exits 0:

```
pool 35,055  scheduled 5,208  proteins 3,798  slots 2,301  rows 2,301  mode Mtm  -> C:\out\cadenza_targets.csv
```

`pool` is precursors passing the q-value cutoff. `scheduled` is targeted
precursors (passengers are not counted). `proteins` is protein groups
covered. `slots` is MTM or PRM slots. `rows` is CSV data rows: one per slot
in MTM, one per precursor in PRM.

| Exit code | Meaning |
|---|---|
| 0 | CSV written |
| 1 | Unexpected failure, for example an unreadable report or bad JSON (`error: <Type>: <message>` on stderr) |
| 2 | Nothing could be scheduled with these settings; no CSV written |
| 64 | Wrong number of arguments |
| 66 | Report or PARAMS file not found |

### Defaults

Unset fields take the **app's** defaults, meaning what a fresh Cadenza
window uses before anyone touches a control. These match the Core
`SchedulingParameters` defaults except for two settings, marked below. The
older `schedule` and `schedule-diann` verbs keep the plain Core defaults.

| JSON field | Default | Notes |
|---|---|---|
| `mode` | `Mtm` | `Mtm` or `Prm` |
| `objective` | `MaximizeProteins` | `Balanced`, `MaximizeProteins`, `MaximizePeptides` |
| `isolationWindowTh` | 3.0 | Max MTM window width (Th) |
| `prmIsolationWidthTh` | 0.7 | PRM / solo-slot width (Th); added to the member span for the CSV isolation width |
| `fragmentTolDa` | 0.5 | |
| `cycleBudget` | 100 | Concurrent slots per cycle |
| `firingPadSec` | 15 | Pad on each side of the peak (s) |
| `rtBinMin` | 0.05 | |
| `qValueCutoff` | 0.01 | DIA-NN precursor q-value |
| `minPeptidesPerProtein` | 1 | |
| `maxPeptidesPerProtein` | 5 | |
| `peptideRankCap` | 0 | 0 = no cap |
| `minTargetIonShare` | **0.015** | App default. Core default is 0 |
| `maxTargetsPerSlot` | 0 | 0 = no cap |
| `backgroundFragmentClash` | **true** | Always on in the app. Core default is false |
| `chargeHandling` | `SameChargePerSlot` | or `AllowMixed` |
| `proteinRanking` | `SummedIntensity` | `SummedIntensity`, `ProteinQValue`, `ProvidedListOrder` |
| `peptideRanking` | `PrecursorQuantity` | or `QValue` |
| `normalizedCollisionEnergy` | 28.0 | Written to the CSV's HCD CE column only |
| `forcedPeptides` | none | Array of entries, same syntax as the app's box (`PEPTIDE`, `PEPTIDE/2`) |
| `targetProteins` | none | Array of accessions or gene names, resolved against the report as the app does |
| `targetMode` | `FirstThenFill` | or `Exclusive` |

Field names are camelCase and enums are strings. Example:

```json
{
  "mode": "Mtm",
  "cycleBudget": 80,
  "normalizedCollisionEnergy": 27,
  "targetProteins": ["P02768", "APOA1"],
  "targetMode": "FirstThenFill"
}
```

## Interactive: GUI launch options

```
SkylineCadenza.exe --report <report.parquet> [--export-dir <folder>]
```

`--name value` and `--name=value` are both accepted.

- If the first argument starts with `--`, Cadenza runs **standalone**. It
  makes no Skyline connection, and **Push targets to Skyline** stays
  disabled. The Skyline status line says why.
- `--report` selects the **Empirical (DIA-NN report)** source, fills in the
  report path, and runs **Load + schedule** as soon as the window opens.
  Every other setting starts at the app's defaults (the table above).
- `--export-dir` makes **Export Thermo CSV...** open in that folder,
  creating it if needed. The default file name is still
  `cadenza_targets.csv`, and the person can rename it. A caller watching the
  folder should pick up any new `*.csv`.
- If an option is unknown or has no value, the error appears in the status
  line, along with the usage.

Launching from Skyline's Tools menu is unchanged. There `args[0]` is the
pipe name from `$(SkylineConnection)`, which never starts with `--`.
Launching with no arguments also behaves as before and attaches to the most
recent running Skyline.

## Tests

`ExportThermoCliTests` writes a small synthetic DIA-NN report, runs
`export-thermo` in both MTM and PRM, and checks that the CSV equals
`ThermoCsvWriter` output from the app's Load pipeline. When the reference
Stellar report is on disk, the test repeats that check on it.
`LaunchOptionsTests` covers argument parsing.
