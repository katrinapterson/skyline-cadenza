# Skyline integration

How Cadenza connects to a running Skyline document, what it reads from the
document, and what it writes back. The user-facing surface is the
"Connect to Skyline" / "Load from Skyline document" / "Update Skyline
document" buttons in the Cadenza UI.

## Transport: JSON-RPC over a named pipe

Cadenza talks to Skyline through Skyline's own JSON-RPC server (the same
transport `SkylineMcpServer` uses upstream). When Skyline launches an
external tool from `Tools > External Tools`, it passes the pipe name as
`args[0]`. The pipe name Skyline hands out is the LEGACY ToolService pipe
name; the JSON server listens on a derived name with the
`SkylineMcpJson-` prefix. `SkylineSession.FromArguments` transforms the
incoming `args[0]` to the JSON pipe name before connecting; connecting to
`args[0]` directly hits the legacy binary server and the first read comes
back as `JsonReaderException: 0x00 is invalid start of value`.

The session is a connect-per-call factory, NOT a long-lived connection.
Every call to `SkylineSession.Execute` opens a fresh pipe, makes the
request, and disposes. This matches Skyline's upstream usage; holding the
pipe open across calls works for the first request and then fails with
`0x00 is invalid start of value` after idle because Skyline expects each
exchange to be self-contained.

If `args[0]` is missing (developer launching the .exe standalone), the
session falls back to scanning `~/.skyline-mcp/connection-*.json` for any
running instance and uses the most recently connected one.

## What Cadenza reads from the document

When the user clicks "Load from Skyline document", `SkylineLibraryLoader`
runs a single custom report against the live document to fetch one row per
`(peptide, charge, fragment)`. From those rows it materialises a
`Candidate` list:

- Peptide identity: `PeptideModifiedSequence`, `PrecursorCharge`,
  `PrecursorMz`, `ProteinName`. Modification syntax is normalised through
  `NormalizeModifiedSequence` so DIA-NN, Carafe, and Skyline forms collapse
  to the same key.
- Fragment list: `ProductMz`, `ProductCharge`, `LibraryIntensity`. The
  top-6 fragments by intensity per `(peptide, charge)` survive into the
  candidate's `Fragments` array; the top-4 by m/z become the
  `Top4Fragments` the scheduler uses for the fragment-clash check.
- Retention time: peptide-level only. Per-replicate columns
  (`PredictedResultRetentionTime`, `BestRetentionTime`) multiply rows by
  the replicate count and blow up on documents with many imports. The
  loader probes `PredictedRetentionTime`, `LibraryRetentionTime`,
  `AverageMeasuredRetentionTime`, `BestRetentionTime`,
  `PeptideRetentionTime`, `RetentionTime` in that order and picks the
  first whose sample of 200 values contains a non-zero entry. The user
  can override the column from the UI if none of the candidates work.
- Peak boundaries: pulled from the active BLIBs Skyline has attached to
  the document. `SkylineBlibDiscovery` reads the document's library list
  (BLIBs by priority order), and `BlibRetentionTimeReader` pulls per-
  peptide `startTime` and `endTime` from each BLIB's `RetentionTimes`
  table. Cadenza takes the union of every replicate's boundary (widest
  start, widest stop) as the candidate's `RtStart` / `RtStop`. The
  scheduler's time padding (a drift buffer) is then layered on top,
  NOT on the peak
  shape itself but on cross-acquisition drift between the source data
  and the new run.

  If no BLIB row exists for a peptide, the boundary is synthesised as
  `RT +/- peakHalfWidthMin` and the candidate is reported in
  `PeptidesFromSynthesizedBoundaries`.

The loader also caches the BLIB paths it consulted on the result; the UI
prints them on the status bar so the user can sanity-check that the
peak-boundary source was what they expected.

## What Cadenza writes back

The "Update Skyline document" button runs `UpdateSkylineAsync`, which
does four things in sequence. Each step is gated on the previous one
succeeding.

### 1. Save guard

Cadenza refuses to push if the Skyline document hasn't been saved. The
assay BLIB needs to live at a stable on-disk path next to the `.sky` file
so Skyline's library registration sticks across saves; writing to `%TEMP%`
made it look like the registration had worked (Skyline reports success on
the `RunCommand` call) but the registration silently dropped the next
time the document was saved. The user is told to save first and retry.

### 2. BLIB write

`BlibAssayWriter.Write` emits a self-contained BiblioSpec
non-redundant `.blib` next to the user's `.sky` file:
`Cadenza-yyyyMMdd-HHmmss.blib` (local time, since assay names are
per-workstation unique and a local timestamp is more readable to the
user). The BLIB is schema version 8 (the version that introduced
`RefSpectra.startTime` / `endTime`):

- `RefSpectra`: one row per scheduled peptide, holding the top-6 library
  fragments as the reference peak list (zlib-compressed BLOBs of
  little-endian doubles for m/z and floats for intensity), plus the
  candidate's RT apex and the per-replicate scheduling window
  (`startTime` / `endTime` from the BLIB or the scheduler's padded
  window if synthesized).
- `Modifications`: UniMod IDs in the modified sequence are translated to
  monoisotopic mass deltas before being written. BiblioSpec's
  `LibKeyModificationMatcher` chokes on `C[UniMod:4]` with `FormatException`
  ("The number 'UniMod:4' is not in the correct format"); writing
  `C[+57.0]` instead avoids the matcher's int-parse path. The mass
  deltas come from a hardcoded UniMod table covering Carbamidomethyl
  (4), Acetyl (1), Phospho (21), Oxidation (35), and several others.
- `Proteins` + `RefSpectraProteins`: the protein-group identifier per
  scheduled peptide, so Skyline's library viewer can roll spectra up to
  proteins.

The `SqliteConnection` uses `Pooling=False`. `Microsoft.Data.Sqlite`
pools connections by default, which keeps the underlying file handle
open on Windows even after the `using` block disposes the connection.
Tests that delete the file straight after the write hit
`IOException: file is being used by another process` without this.

### 3. Library registration

The BLIB is attached to the running document via a single `RunCommand`
call:

```
--add-library-path=<blib path>
--add-library-name=<assay name>
```

This routes through Skyline's `CommandLine.SetLibrary(name, path)`, which
calls `LibrarySpec.CreateFromPath` to dispatch on the file extension
(`.blib` -> `BiblioSpecLiteSpec`) and appends to
`PeptideSettings.Libraries.LibrarySpecs`. Cadenza initially tried the
JSON-RPC `AddSettingsListItem("Libraries", ...)` path, which fails because
Skyline's `JsonToolServer.DeserializeSettingsItem` looks for a static
`Deserialize(XmlReader)` on the generic argument of `SettingsListBase<T>`,
and for the Spectral Libraries list `T` is the abstract `LibrarySpec`
(which has no such method, only its concrete subclasses do). The
`RunCommand --add-library-*` path is the routes-around-it workaround.

The `--add-library-path` is sent BEFORE `--add-library-name` in the same
batch. They're documented as order-insensitive (SkylineCmd parses all args
before running handlers), but in practice the library was sometimes
leaving the active list unattached when name came first. Path-first turns
out to be stable.

### 4. Settings configurator

`SkylineSettingsConfigurator.ConfigureAsync` derives a recommendation
from the scheduled candidates (precursor charge set, fragment charge set,
peptide length range, ion types y/b/p, MS1 isotope filter, RT filter
tolerance) and applies it through `RunCommand` calls. Each flag goes in
its own batch except for the two mutual-validation cases:

- DIA acquisition method + `Results only` isolation scheme are bundled
  because `TransitionFullScan.DoValidate` rejects each one alone:
  isolation-scheme without DIA throws "No other full-scan MS/MS filter
  settings are allowed when precursor filter is none"; DIA without
  isolation-scheme throws "An isolation window width value is required
  in DIA mode".
- MS1 precursor isotope filter (`Count` + threshold + analyzer + ppm)
  goes in one batch because Skyline runs `ValidateRes` against the
  combined end state.

The configurator's `Recommendation` record is also stamped into the
status bar so the user can see exactly which settings Cadenza pushed.
The `ToUiInstructions` fallback writes the same recommendation as
click-by-click Skyline UI instructions in case a flag fails on an
older Skyline version.

### 5. Transition-list import

`PeptideTransitionListBuilder.Build` emits a peptide-style transition
list CSV (header: `Protein Name`, `Peptide Modified Sequence`,
`Precursor Charge`, `Precursor m/z`, `Product m/z`, `Product Charge`,
`Explicit Retention Time`, `Explicit Retention Time Window`, `Note`).
Each scheduled precursor produces:

- 3 precursor isotope rows (M+0, M+1, M+2). Skyline recognises a
  transition where `Product Charge == Precursor Charge` and the product
  m/z is at the precursor's monoisotope spacing as a precursor
  transition; mass step is the 12C -> 13C neutron offset
  (1.003355 Da / charge per isotope step).
- Up to 6 fragment rows (`FragmentsPerPeptide == BlibAssayWriter.PeaksPerSpectrum`),
  each carrying its real `FragmentIon.Charge`. An earlier build
  hardcoded `Product Charge = 1` and Skyline rejected the +2 y/b
  fragments DIA-NN selects on longer tryptic peptides with
  "no matching product ion".

Heavy isotope-labeled precursors (e.g. spiked PRTC standards) are
**omitted** from the transition list. Skyline reports a heavy precursor's
m/z alongside the bare light sequence (the isotope label is not part of
`PeptideModifiedSequence`), so a naive row pairs a light sequence with a
heavy m/z. Skyline's headless transition-list importer can't apply the
document's `explicit_decl` heavy modification ("isotope labeled entry
found without matching settings") and, because it uses an early row to
validate the precursor-m/z column, a single such row aborts the whole
import with the misleading "No valid precursor m/z column found, line 1".
`SkylineLibraryLoader` captures each precursor's `IsotopeLabelType`, and
`PeptideTransitionListBuilder` skips any precursor that isn't light. Those
peptides still reach the assay as heavy through the registered spectral
library (step 3). The number omitted is reported on the status line.

The CSV is written to `%TEMP%` and imported with
`RunCommand --import-transition-list=<path> --ignore-transition-errors`.
`--ignore-transition-errors` degrades any remaining per-row rejection to a
warning rather than aborting the import, so one unexpected row can't
discard the entire scheduled assay. The import is slow
(~10 transition rows/sec) because Skyline's `MassListImporter` rebuilds
the document tree per protein; this is on the Skyline team to address.
Cadenza shows a heartbeat status every 10 s with the elapsed time and
ETA so the UI doesn't look hung. The temp CSV is deleted in `finally`.

## Failure modes the user might see

- "Save the Skyline document first ...": step 1 guard.
- "BLIB written ... but registering it with Skyline failed: ...": step 3
  reported a non-empty stdout containing "Error" or "conflict". The
  message includes Skyline's verbatim output and tells the user to add
  the library manually via Peptide Settings then re-run the push.
- "library registration failed: RunCommand --add-library-* threw: ...":
  the JSON-RPC pipe itself errored. Usually means Skyline crashed or the
  pipe closed; reconnect from the Skyline button and try again.
- Library appears registered but doesn't show up in the Peptide Settings
  dialog: known cosmetic bug in Skyline's settings dialog cache. The
  library IS attached (the document XML and a fresh "Update Skyline
  document" push both confirm it); closing and reopening the Peptide
  Settings dialog forces a refresh. Nick Shulman has a fix planned for
  a future Skyline release.
