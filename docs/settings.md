# Settings reference

Every user-tunable knob exposed by the Cadenza UI maps to a property on
the `SchedulingParameters` record (`src/SkylineCadenza.Core/Scheduling/
SchedulingParameters.cs`). Each property has a one-to-one
`[ObservableProperty]` on `MainViewModel`, a XAML control, and a
`partial void On...Changed` hook that re-schedules whenever the user
moves a knob.

This page documents what each setting controls, its default value, and
which other settings it interacts with. For the algorithm itself see
[`coverage-objectives.md`](coverage-objectives.md); for what gets pushed
back to Skyline see [`skyline-integration.md`](skyline-integration.md).

## Acquisition mode

### Mode (PRM | MTM)

- Default: `MTM`.
- PRM = one precursor per slot; the quadrupole window is
  `PrmIsolationWidthTh` wide on every fired precursor.
- MTM = multiplexed; multiple precursors per slot, slot width is
  `IsolationWindowTh`.
- Cadenza pushes the matching `--full-scan-acquisition-method` flag
  (PRM or DIA) into Skyline; MTM additionally sets isolation scheme
  to `Results only`.

### Charge handling (`SameChargePerSlot` | `AllowMixed`)

- Default: `SameChargePerSlot`.
- `SameChargePerSlot`: MTM slots may only contain precursors of
  identical charge. Cleaner collision-energy assignment and simpler
  downstream interpretation.
- `AllowMixed`: MTM slots may contain mixed charges. The Thermo CSV
  writer reports each slot's majority charge (ties resolve to the lower
  z, since +2 is the canonical tryptic default).
- Ignored in PRM mode (each slot has one precursor).

### Normalized Collision Energy (NCE)

- Default: 28.0.
- Emitted into the Thermo inclusion CSV, NOT used by the scheduler.

## Slot geometry

### Isolation window width (Th)

- Default: 3.0 Th.
- The nominal slot width in MTM mode. Each MTM slot's center span
  (max member m/z minus min member m/z) is capped at
  `IsolationWindowTh - PrmIsolationWidthTh` so every member's
  quadrupole window fits fully inside the slot's nominal isolation
  window (the slot-edge rule).
- Ignored in PRM mode.

### PRM quadrupole isolation width (Th)

- Default: 0.7 Th.
- Used in three places:
  1. PRM mode: the isolation width for every scheduled precursor.
  2. MTM solo slots (slots that ended up with one member): the
     instrument isolation written into the Thermo CSV.
  3. MTM slot-edge rule (above): every member's center must sit at
     least `PrmIsolationWidthTh / 2` from each slot edge.

### Fragment tolerance (Da)

- Default: 0.5 Da.
- Used by the fragment-clash check when joining a peptide to an existing
  slot. Two slot members' top-4 fragments are "clashing" if any pair is
  within this tolerance. Tight tolerance = more aggressive clash
  rejection = more singleton slots = lower multiplexing efficiency.

## RT budget

### Cycle budget (concurrent slots per cycle)

- Default: 100.
- Maximum concurrent acquisition slots at any single RT bin. This is the
  hard load-balancing constraint: the scheduler will refuse to add a
  peptide whose padded scheduling window would push any of its RT bins
  past this number. (The scheduling window is the peak boundary plus
  the time padding above.)
- The budget is per RT bin, not global. A protein with peptides spread
  across the gradient at distinct RTs consumes one slot at each of
  those bins only.

### Time padding (sec)

- Default: 15.0 s (0.25 min).
- Extra time added on each side of every peak's RT range when the
  scheduler builds the scheduling window. This is a drift buffer:
  the source data's peak boundaries already include per-replicate
  peak-shape variance, so this padding is for retention-time drift
  between the source acquisitions and the new run, not for the peak
  shape itself.
- The matching code symbol is `SchedulingParameters.FiringPadSec`
  (kept unrenamed for now to avoid churning the public Core API; the
  UI label and docs use "time padding").

### Predicted peak half-width (sec, Carafe workflow only)

- Default: 6.0 s (0.10 min).
- Half-width of the synthesised peak on the Carafe path, where the
  predicted library carries only an apex (`Tr_recalibrated`). The
  candidate's unpadded peak is `apex +/- this`, and the time-padding
  buffer above is added on top. Empirically (DIA-NN search of the
  matched run) real peaks are wider - median full width ~0.40 min, i.e.
  ~0.20 min (12 s) half-width - so the default is on the tight side;
  raise it toward 10-12 s for more realistic windows.
- Ignored on the DIA-NN and Skyline-document workflows, which use the
  real measured `RT.Start` / `RT.Stop` (or BLIB) boundaries.
- The matching code symbol is the `peakHalfWidthMin` parameter of
  `CandidateBuilder.BuildFromCarafe`; the UI exposes it in seconds.
  Adjusting it rebuilds the Carafe candidate pool and reschedules live.

### RT bin width (min)

- Default: 0.05 min (3 s).
- Width of the RT occupancy bins used for the per-RT-bin budget check.
  Smaller bins = finer-grained load curve = slightly higher schedule
  density but more memory and CPU. Larger bins = coarser smoothing and
  spurious budget rejections near peak edges.

## Q-value filter (DIA-NN ingest)

### Q-value cutoff

- Default: 0.01.
- DIA-NN `Q.Value` cutoff applied at the candidate-build stage. Rows
  above this cutoff are dropped before the scheduler sees them.

## Coverage objective

### Objective (`Balanced` | `MaximizeProteins` | `MaximizePeptides`)

- Default: `MaximizeProteins`.
- Selects the cover-pass strategy AND the load-up cap. See
  [`coverage-objectives.md`](coverage-objectives.md) for details.

### Min peptides per protein

- Default: 1.
- Hard filter, applied AFTER the load-up loop finishes. Groups whose
  final scheduled count is below Min are dropped from the assay
  entirely (their slots are freed and the load curve is recomputed).
  Setting Min > 1 means "every protein must reach this depth or get
  cut from the assay"; proteins that can't reach it disappear.
- Min does NOT participate in the load-up cap. The load-up loop
  always uses Max as the per-protein ceiling (under Balanced and
  MaximizeProteins) regardless of Min.

### Max peptides per protein

- Default: 5.
- Per-protein ceiling for the load-up loop. The load-up loop is
  breadth-first round-robin: lap 1 brings every covered protein up
  to 2 peptides where possible, lap 2 to 3, and so on. A protein
  drops out of the loop on the lap where it either hits Max or runs
  out of candidates that fit, whichever comes first. Proteins with
  only 1 or 2 viable peptides stop early and don't block the others
  from continuing to grow.
- Honoured by both `Balanced` and `MaximizeProteins`. Under
  MaximizeProteins, once every coverable protein has its first
  peptide the load-up loop fills toward Max using whatever cycle
  budget remains.
- Ignored under `MaximizePeptides`: the load-up loop keeps going,
  lap by lap, until either the cycle budget is saturated or every
  protein has exhausted its candidate queue.
- Cover-pass behaviour is independent of Max under all objectives, so
  changing Max does not change the protein-coverage curve, only the
  per-protein peptide depth.

## Peptide / protein ranking

### Protein ranking (`SummedIntensity` | `ProteinQValue` | `ProvidedListOrder`)

- Default: `SummedIntensity`.
- `SummedIntensity`: sum of DIA-NN `Precursor.Quantity` per group.
  Highest first. Sensible when you want to schedule the most-intense
  proteins first.
- `ProteinQValue`: best (smallest) `PG.Q.Value` per group. Lowest
  first. Sensible when the source report is borderline-noisy and you
  want to prioritise high-confidence ID groups.
- `ProvidedListOrder`: keep the user's target list in insertion order.
  Sensible when the user has a hand-curated priority list.

### Peptide ranking (`PrecursorQuantity` | `QValue`)

- Default: `PrecursorQuantity`.
- `PrecursorQuantity`: DIA-NN `Precursor.Quantity`, descending (highest
  first). This is the published behaviour and the default.
- `QValue`: DIA-NN `Q.Value`, ascending (lowest first). Useful when
  intensity rankings are unreliable (e.g. matrix interference suspected
  to inflate certain precursors).
- Unique-before-razor tier always wins over either of these: a unique
  peptide will always be tried before a razor peptide regardless of
  quantity or q-value.

## Target list

### Target proteins (multi-line text box)

- Default: empty.
- A user-supplied set of protein accessions to focus the assay on. One
  accession per line; the count is shown live next to the box as a
  parse check.
- Empty = no filter; the scheduler considers every protein the report
  sees.

### Target mode (`Exclusive` | `FirstThenFill`)

- Default: `FirstThenFill`.
- `Exclusive`: only proteins in the target list are eligible.
  Off-list proteins are dropped before the scheduler sees them.
- `FirstThenFill`: target-list members are scheduled first (in the
  order they were supplied), then the remaining budget fills with
  off-list proteins in the regular priority order.
- Ignored when the target list is empty.

## Skyline integration (read-side)

### RT column override (text box)

- Default: empty (auto-probe).
- When loading from a running Skyline document, the candidate column
  for peptide-level retention time. Empty = Cadenza probes
  `PredictedRetentionTime`, `LibraryRetentionTime`,
  `AverageMeasuredRetentionTime`, `BestRetentionTime`,
  `PeptideRetentionTime`, `RetentionTime` in order and uses the first
  one with non-zero values in the sample. Set this when the document
  has none of those configured (or all return zero) and the user knows
  which column actually holds the RT.
- See [`skyline-integration.md`](skyline-integration.md) for details.

## How settings get applied

Every property has a `partial void On<Name>Changed` hook on
`MainViewModel` that calls `RescheduleAsync()`. Moving any knob
re-runs the scheduler on the current candidate pool and re-renders
every plot. The "Update Skyline document" button does NOT re-run the
scheduler; it pushes whatever is currently displayed.

When adding a new setting:

1. Add the field + XML doc to `SchedulingParameters`.
2. Add an `[ObservableProperty]` backing field to `MainViewModel`.
3. Bind a XAML control to the new property in `MainWindow.xaml`.
4. Add the `partial void On<Name>Changed(...) => _ = RescheduleAsync()`
   hook.
5. Surface the value in `CurrentParameters()` so the scheduler sees it.
6. Cover the behaviour with a test in `SkylineCadenza.Tests`.
