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

### Background fragment clash (always on, no control)

- MTM only. The clash check above only compares a slot's own members. This
  one also covers every identified precursor the quadrupole transmits with
  them, including never-scheduled background. It looks at the identified
  precursors that co-elute with a target inside the emitted band
  (`MzMin - PrmIsolationWidthTh / 2` to `MzMax + PrmIsolationWidthTh / 2`)
  and have a top-4 fragment within the fragment tolerance of one of the
  target's top-4 fragments. Their signal on that fragment is added up, each
  weighted by how much of its peak co-elutes with the target. The fragment
  counts as interfered once that total reaches 1% of the fragment's own
  signal (a fixed tolerance, not a setting), so a background precursor 1000x
  weaker than the target does not block it. All four quantified fragments
  must be clean, and a refused target falls through to its protein's next
  peptide.
- Signal per fragment is `PrecursorQuantity * Intensity / sum(Intensity)`.
  Only identified precursors (quantity > 0) are known, so this is a lower
  bound on the real interference.
- Measured on the 08/03 HeLa GPF report (App defaults, truncated-FASTA
  artefacts excluded; expected re-identification from a model fitted on the
  acquired MTM run). MaximizeProteins went from 5,777 to 5,745 targets, from
  ~5,090 to ~5,109 expected re-identified targets (88.1% to 88.9%), and from
  ~3,865 to ~3,851 expected recovered proteins. Balanced went from 5,013 to
  5,039 targets and from ~4,547 to ~4,592 expected re-identified (90.7% to
  91.1%). Targets with a background hit on a top-4 fragment fell from about
  24-28% of the schedule to 13-16%.
- Code symbol: `SchedulingParameters.BackgroundFragmentClash` (false in
  Core so the golden fixtures are unchanged; the App always sets it true).

### Max targets per slot

- Default: 0 (no cap).
- Maximum number of targeted precursors (drivers) one MTM slot may hold. A
  join into a full slot is refused, so the target goes to another slot,
  opens a solo slot, or falls through to its protein's next peptide.
- Untargeted co-isolated precursors and passengers do not count.
- Enforced alongside the co-isolation gates (join, look-ahead join check,
  shrink-and-place; a new solo slot always holds one target).
- Code symbol: `SchedulingParameters.MaxTargetsPerSlot`.

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

### Carafe detectability predictions + cutoff (Carafe path only)

- **Detectability predictions**: one or more Carafe `*_detectability_pred.tsv`
  files (multi-select) or a folder, entered separately from the Bromo
  peptide-ranking overlay. Parsed by `PredictedRankingReader`. Two roles:
  1. the cutoff filter below, and
  2. the within-protein ranking **fallback** - when no Bromo file is loaded,
     detectability drives peptide priority (`PeptidePriority.PredictedRanking`)
     instead of the no-signal `PrecursorQuantity` stand-in. Bromo, if supplied,
     takes precedence for ranking; detectability then only filters.
- **Detectability cutoff** (default **0.5**): precursors whose looked-up Carafe
  detectability is below this are dropped from the candidate pool before
  scheduling. Precursors with no prediction are kept (unknown is not filtered).
  Applied only when a detectability file is supplied. Empirically ~0.5 removes
  peptides that are ~90% unlikely to be observed while retaining ~87% of real
  detections; raise toward 0.9 for high-confidence-only targets (at the cost of
  recall). Wired to `CandidateBuilder.BuildFromCarafe(..., detectability,
  minDetectability)`; moving the slider rebuilds the Carafe pool and reschedules
  live.

### Predicted peak half-width (sec, Carafe path only)

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

## DIA-NN + Carafe backfill

### Backfill unused schedule time with Carafe predictions (checkbox)

- Default: **off**. Requires both a DIA-NN report and a Carafe library.
- **On**: builds a two-tier candidate pool. The confirmed DIA-NN targets
  (experimental tier) are scheduled to exhaustion first; the Carafe-predicted
  targets (predicted tier) then reuse whatever cycle budget is left. Predicted
  precursors that duplicate a confirmed one on (stripped peptide, charge) are
  dropped. Predictions are scoped to the accessions DIA-NN observed plus any
  target-protein list, so they add depth to confirmed proteins and cover
  proteins DIA-NN missed.
- **Off**: the DIA-NN path is unchanged - the Carafe library only overrides the
  measured fragments and adds no precursors.
- The predicted tier honors the same **cycle budget** and **max peptides per
  protein** as the confirmed tier (a protein already at its Max from
  experimental peptides gets no predicted peptides), and its breadth-vs-depth
  follows the selected **coverage objective**. The **peak half-width**,
  **detectability**, and **Bromo ranking** controls apply to the predicted tier;
  the status line reports the confirmed + predicted split.
- Because predicted precursors join the pool, **protein ranking** (summed
  intensity) and parsimony grouping see the combined evidence, so the confirmed
  schedule can shift slightly versus DIA-NN alone. Confirmed coverage is not
  reduced. Takes effect on the next Load.
- Code path: `MainViewModel.LoadCombinedBackfillAsync` ->
  `CandidateBuilder.BuildCombined`; the two-tier scheduling lives in
  `Scheduler.Run` keyed on `Candidate.Source` (`CandidateSource`).

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

## Forced peptides

### Forced peptides (multi-line text box + "Open peptide list...")

- Default: empty (no effect).
- Peptides placed before anything else, so the same validation or
  intentionally targeted peptides can be carried across different
  schedules. The normal cover and load-up passes then fill the rest of the
  schedule around them.
- One entry per line: `PEPTIDEK` (any charge), `PEPTIDEK+2`, `PEPTIDEK/2`,
  `PEPTIDEK++`, or a charge in a second comma / tab / space separated column.
  Modified sequences in Skyline, DIA-NN or Carafe syntax must match the
  modifications (masses compared to two decimals); an unmodified entry
  matches every modified form. Blank lines, `#` comments and header rows are
  skipped.
- Forced peptides obey every quality and instrument rule: the ion-share
  floor, the background fragment clash, the member fragment clash, window
  geometry and the cycle budget. They ignore the breadth rules: ranking, the
  rank cap, max and min peptides per protein, and an exclusive target list.
- Each forced precursor is placed the normal way under the schedule's own
  parameters: join a compatible window, else open a new one, else
  shrink-and-place. So a core set carried across schedules is scheduled the
  way each schedule's settings dictate, which is what makes the schedules
  comparable. A peptide without a charge tries its best-ranked charge state
  first. If that fails the rules, it falls back to its other charge states,
  in rounds, so no entry's fallback takes budget another entry's best charge
  needed. The fill then treats forced windows like any other window: other
  targets may join or reshape them under the same rules, but the fill's
  shrink-and-place never drops a forced precursor.
- The status line reports how many entries were scheduled (and how many on
  a fallback charge), which are not in the library (no identification, so
  no RT), and which failed the rules.
- Measured on the 08/26 HeLa GPF report: 30 random peptides the normal
  schedule did not pick were forced. 28 were scheduled (1 on a fallback
  charge). The other 2 failed the background fragment clash in their own
  solo windows. The fill still reached 5,093 targets and 3,704 proteins,
  against 4,915 and 3,688 without forcing.
- Code symbols: `SchedulingParameters.ForcedPeptides`,
  `ForcedPeptides.Parse`, `ScheduleResult.ForcedOutcomes`.

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
