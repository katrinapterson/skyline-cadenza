# Ingest sources

Cadenza accepts a peptide candidate pool from three different sources.
Each source contributes a different subset of the fields a `Candidate`
needs, so the path the data takes through the ingest layer is different
for each. This page documents what each source supplies, where the gaps
get filled, and how to interpret the result when a field is synthesised
or proxied.

## What a Candidate needs

The scheduler reads these fields off each `Candidate`:

- Identity: `PrecursorId`, `StrippedSequence`, `ModifiedSequence`,
  `PrecursorCharge`, `PrecursorMz`.
- Peak: `RtApex`, `RtStart`, `RtStop` (unpadded; the scheduler adds
  its time padding on top as a drift buffer).
- Quality: `QValue`, `ProteinQValue`. Used for filtering and for
  protein-ranking when `ProteinPriority.ProteinQValue` is selected.
- Quantity: `PrecursorQuantity`. Used for peptide-ranking and for the
  `ProteinPriority.SummedIntensity` protein order.
- Parsimony: `ProteinGroup`, `PeptideType` (`unique` | `razor` | other).
  Drives the smallest-groups-first cover order and the
  unique-before-razor tie-break.
- Fragments: `Fragments[]` (top-N by intensity) and the derived
  `Top4Fragments` (top-4 by m/z) for the slot-join fragment-clash
  check.

Anything not in the source has to be synthesised from elsewhere or
defaulted. The way each source fills those slots is described below.

## 1. DIA-NN report

Path: open `report.parquet` from a DIA-NN run.

`DiannParquetReader` streams rows for every (peptide, charge, run) and
materialises:

- Identity from `Modified.Sequence` + `Precursor.Charge`. Stripped
  sequence is derived by removing the bracket/paren UniMod annotations.
- Peak: `RtApex` from `RT`; `RtStart` / `RtStop` from DIA-NN's
  `RT.Start` / `RT.Stop` (the per-replicate peak boundaries DIA-NN
  records in the report).
- Quality: `Q.Value` (peptide-precursor q-value) and `PG.Q.Value`
  (protein-group q-value).
- Quantity: `Precursor.Quantity` (DIA-NN's real measured intensity).
- Parsimony: not in DIA-NN. Cadenza runs its own parsimony pass on
  the (peptide, protein) edges DIA-NN reports in `Protein.Group`,
  building unique/razor assignments per peptide.
- Fragments: parsed from DIA-NN's `Fr.N.Id` columns ("y6^2/447.71" =
  y6, +2, m/z 447.71). Charge defaults to +1 when the source row's
  `^z` field is absent. These are REAL measured fragments, not
  predictions.

A row is dropped at this stage if its `Q.Value` exceeds the user's
`QValueCutoff` (default 0.01) or if its peptide isn't assigned to any
parsimonious group.

DIA-NN ingest is the most data-rich path: every field comes from a real
measurement.

## 2. Carafe AI-predicted spectral library

Path: open `carafe-library.tsv` (the AI-predicted DIA-NN-format
spectral library Carafe produces from a target protein list).

`CarafeLibraryReader` streams the TSV, early-rejecting on the
`ProteinID` column, and materialises one `CarafePrecursor` per
(peptide, charge):

- Identity from `ModifiedPeptide` / `StrippedPeptide` /
  `PrecursorCharge`. Modification syntax is Carafe's
  `_C[UniMod:4]DIVIEK_` form; `NormalizeModifiedSequence` collapses it
  to Skyline's `C[UniMod:4]DIVIEK` form downstream.
- Peak: `RtApex` from `Tr_recalibrated` (Carafe's predicted apex).
  `RtStart` / `RtStop` are SYNTHESISED in `CandidateBuilder.BuildFromCarafe`
  as `RtApex +/- peakHalfWidthMin` (default 0.10 min / 6 s), exposed as the
  "Predicted peak half-width" UI slider; adjusting it rebuilds the Carafe
  pool and reschedules live. The scheduler then layers its time padding on
  top. Empirical peaks (DIA-NN search of the matched run) are wider - median
  ~0.40 min full width - so the default is on the tight side.
- Quality: `QValue = 0.0` (there's no q-value on a prediction; the
  filter is upstream, in Carafe's own training data). `ProteinQValue =
  NaN`.
- Quantity: NOT PROVIDED. Carafe predicts fragment relative intensity
  for each (peptide, charge), not a precursor or peptide abundance.
  See below for what Cadenza does with the
  `PrecursorQuantity`-ranked path in Carafe-only mode and how a
  predicted-ranking overlay supplies a real within-protein order.
- Parsimony: not in Carafe. Cadenza builds a peptide-to-proteins map
  from the Carafe `ProteinID` column (semicolon-separated
  UniProt-style `sp|ACC|NAME` ids are parsed to bare accessions) and
  runs its own parsimony pass.
- Fragments: from the per-row `FragmentMz` / `RelativeIntensity` /
  `FragmentCharge` columns (when `FragmentCharge` is absent, charge
  defaults to +1). These are PREDICTED fragments with predicted
  relative intensities (normalised within a single (peptide, charge)
  spectrum). Cadenza keeps up to 12 per (peptide, charge),
  intensity-sorted descending, then takes the top-4 by m/z for the
  slot-join clash check.

### Peptide ranking in Carafe-only mode

Carafe gives you predicted fragment relative intensities per
(peptide, charge), normalised within each spectrum (the fragments of
one peptide sum to a fixed constant: 100, or 1, depending on
convention). There is no precursor abundance, no peptide-level
quantity, and no across-peptide intensity calibration.

When no predicted-ranking overlay is supplied, the code stuffs the
sum of the top-N fragment relative intensities into the candidate's
`PrecursorQuantity` field so that `PeptidePriority.PrecursorQuantity`
has SOMETHING to sort on, but that sum carries no useful ranking
information:

- Within a single peptide, the relative intensities sum to the
  normalisation constant by construction. The top-N sum is
  bounded above by that constant, so it can only differ between
  peptides by how concentrated the spectrum is in its top N
  fragments. Concentration is uncorrelated with sensitivity or
  detectability in any directly useful way; if anything it favours
  peptides whose predicted spectra are sparse (a peptide with
  exactly N or fewer fragments saturates the top-N sum at the
  normalisation constant).
- Between peptides, there's no quantitative signal at all because
  the per-spectrum normalisation wipes out the absolute scale.

So without an overlay, peptide ranking on the Carafe path falls back
on the unique-before-razor tier (which IS meaningful), and the
`PrecursorQuantity` tie-break inside each tier is effectively noise.

#### Predicted-ranking overlay (recommended)

To get a real WITHIN-protein peptide ranking, supply a ranking source
alongside the Carafe library. There are two inputs, both read by
`PredictedRankingReader` and matched onto the candidate pool by
(peptide, charge) in `CandidateBuilder.BuildFromCarafe`, populating
`Candidate.PredictedPriority`. When either supplies a priority the
peptide-ranking selector auto-switches to
`PeptidePriority.PredictedRanking`, which orders peptides within each
protein by that score, descending (candidates with no match sort last).

The two inputs and their precedence:

- **Bromo peptide ranking** (UI: "Bromo peptide ranking";
  `MainViewModel.PredictedRankingPath`) - a single
  `peptide_rankings.tsv` (`protein, peptide (SEQ|charge),
  mean_pred_score, wins, rank`; score = `mean_pred_score`, higher =
  better, matched on stripped sequence + charge). When present it is the
  ranking source.
- **Carafe detectability predictions** (UI: "Carafe detectability
  predictions"; `MainViewModel.DetectabilityPath`) - one or more
  `*_detectability_pred.tsv` files or a folder
  (`sequence, charge, ModifiedPeptide, detectability_pred, ...`; score =
  `detectability_pred`, matched on Carafe-form `ModifiedPeptide` + charge
  with stripped fallback; shards are merged, higher score wins). This
  input drives the detectability cutoff filter, and **also serves as the
  ranking fallback**: when no Bromo file is loaded, detectability
  populates `PredictedPriority` so peptides are ordered by detectability
  rather than the no-signal `PrecursorQuantity` stand-in.

So the ranking-priority source is: Bromo if provided, else Carafe
detectability if provided, else the `PrecursorQuantity` stand-in.

This is a within-protein signal only. It does NOT give cross-protein
abundance, so it does not feed protein ranking (see below).

What we will never get from Carafe alone: absolute or relative
abundance at the PROTEIN level (between proteins in the same sample).
There's no measurement to predict from.

So when running in Carafe-only mode TODAY:

- Peptide ranking with `PredictedRanking` (auto-selected when an
  overlay is loaded) orders peptides within each protein by the
  predicted detectability / ranking score. This is the recommended
  path. Without an overlay, `PrecursorQuantity` (the default) selects
  on the unique-before-razor tier and then orders randomly within each
  tier, and `QValue` ranking is no better (every candidate has
  `QValue = 0.0`).
- Protein ranking with `SummedIntensity` aggregates that same
  no-signal stand-in across a protein's peptides. With per-spectrum
  normalised intensities, the sum is approximately proportional to
  the protein's peptide count, so it ends up ranking proteins by
  "how many tryptic peptides the protein contributes" rather than
  by abundance. Don't rely on it for protein priority.
- Protein ranking with `ProteinQValue` is meaningless (Carafe's
  `ProteinQValue` is always NaN, which sorts to the end). Don't use
  it.
- Protein ranking with `ProvidedListOrder` is the recommended path
  for Carafe-only mode. Supply a target protein list in the order
  you want the proteins prioritised, and the scheduler will respect
  that order. Most users running Carafe-only have a specific list in
  mind anyway (otherwise they wouldn't be designing a targeted
  assay) so this aligns the tool with the intent.

## 3. Running Skyline document

Path: click "Load from Skyline document" while connected to a running
Skyline document.

See [`skyline-integration.md`](skyline-integration.md) for the
transport details; this section focuses on what the candidate pool
looks like.

`SkylineLibraryLoader` runs one custom report against the live
document to fetch one row per (peptide, charge, fragment), regroups
in memory, and materialises:

- Identity from the document tree.
- Peak boundaries: from the active BLIBs the document has attached.
  `SkylineBlibDiscovery` walks Skyline's library list in priority
  order; `BlibRetentionTimeReader` pulls per-peptide `startTime` /
  `endTime` from each BLIB's `RetentionTimes` table. Cadenza takes
  the union of every replicate's boundary (widest start, widest
  stop). If the peptide isn't in any active BLIB, the boundary is
  synthesised as `RT +/- peakHalfWidthMin` (default 0.30 min on
  this path).
- RT apex: from the document. The loader probes
  `PredictedRetentionTime`, `LibraryRetentionTime`,
  `AverageMeasuredRetentionTime`, `BestRetentionTime`,
  `PeptideRetentionTime`, `RetentionTime` (in that order) and picks
  the first whose sample has a non-zero entry. The user can override
  the column via the UI's "RT column override" box.
- Quality / quantity: not on this path. `QValue` / `ProteinQValue`
  default to 0.0 / NaN; `PrecursorQuantity` defaults to a constant.
  Peptide ranking on a Skyline-document candidate pool falls back to
  the unique-before-razor tier ordering and the alphabetical
  fallback.
- Parsimony: Cadenza runs its own parsimony pass on the
  (peptide, protein) edges the document reports.
- Fragments: from the document. `LibraryIntensity` provides the
  relative intensities that the top-6 selection ranks on;
  `ProductMz` + `ProductCharge` carry the actual fragments.

The Skyline-document path is the "use the document I already have"
case. The candidate pool is shaped by Skyline's own library
configuration; if the document's BLIBs don't carry peak boundaries
for a peptide (e.g. a freshly-built library with no per-replicate
retention info), Cadenza synthesises a tight window and the user
sees the synthesised count on the status bar.

## Combining sources

Two pairings are supported:

- DIA-NN report + Carafe library: DIA-NN provides the candidate
  pool (with its measured RTs, intensities, and q-values), and the
  Carafe fragments override DIA-NN's per-row fragments when a
  matching `(peptide, charge)` exists in the library. This combines
  the strongest signal on each axis: real measured RT and quantity,
  predicted high-quality fragment shortlist.
- Skyline document + Skyline BLIBs: implicit. The "Load from Skyline
  document" path always reads from the document's attached BLIBs.

Carafe + Skyline document is not currently a supported pairing.

DIA-NN + Skyline document is also not a pairing: the two are
alternative entry points to the candidate pool. If you want
DIA-NN's measured data the load-from-Skyline path will not pick it
up; open the parquet directly.
