# Skyline Cadenza Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on `main`. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and `Version` in `src/SkylineCadenza.App/tool-inf/info.properties` is updated to match.

## New Features

- Added a **DIA-NN + Carafe backfill** mode (checkbox "Backfill unused schedule time with Carafe predictions", requires both a DIA-NN report and a Carafe library). Confirmed DIA-NN targets are scheduled first (experimental tier); Carafe-predicted targets then fill whatever schedule capacity is left over (predicted tier), deduped against the confirmed targets on peptide + charge. Predictions are scoped to the proteins DIA-NN observed plus any target-protein list, so they add depth to confirmed proteins and cover proteins DIA-NN missed. Backfill breadth-vs-depth follows the selected coverage objective, and predicted peptides respect the same cycle budget and per-protein cap as confirmed ones (a protein already at its Max from experimental peptides gets no predicted ones). The peak-half-width and detectability controls apply to the predicted tier; the status line reports the confirmed + predicted split. Off (default) leaves the previous behavior, where the Carafe library only overrides DIA-NN fragments and adds no precursors. Note: because predicted precursors join the pool, protein ranking (summed intensity) and parsimony grouping see the combined evidence, so the confirmed schedule can shift slightly versus DIA-NN alone; confirmed coverage is not reduced.
- Carafe detectability now also serves as the within-protein **ranking fallback**: when no Bromo ranking file is loaded, peptide priority is driven by Carafe detectability (auto-selecting `PredictedRanking`) instead of the no-signal fragment-intensity stand-in. Bromo, when supplied, still takes precedence for ranking and detectability then only filters.
- The **Bromo peptide ranking** input is now a single-file, Bromo-only selector (Carafe detectability predictions have their own input).
- Added a separate **Carafe detectability predictions** input (multi-file / folder) plus a **detectability cutoff** slider (default 0.5) for the Carafe path. This is independent of the predicted peptide-ranking overlay, so a user can rank with Bromo and filter by Carafe detectability simultaneously. Precursors below the cutoff are dropped before scheduling (peptides with no prediction are kept); the slider rebuilds the candidate pool and reschedules live. Empirical calibration: ~0.5 removes peptides ~90% unlikely to be observed while keeping most real detections.

- Added a configurable **Predicted peak half-width** slider (seconds) for the Carafe path, replacing the hardcoded 0.10 min synthesised peak. Adjusting it rebuilds the Carafe candidate pool and reschedules live; ignored on the DIA-NN / Skyline paths, which use measured peak boundaries.
- Added a predicted peptide-ranking overlay for Carafe-only mode. Supply a Bromo `peptide_rankings.tsv` or Carafe `*_detectability_pred.tsv` (multi-select several shard files, a folder of shards, or one path per line) via the new "Predicted peptide ranking" input. The score is matched onto the Carafe candidate pool by (peptide, charge) and drives a real within-protein peptide priority (`PeptidePriority.PredictedRanking`), which auto-selects when an overlay is loaded. Replaces the no-signal fragment-intensity-sum stand-in documented in `docs/ingest-sources.md`. Within-protein only; it does not affect cross-protein protein ranking.

## Bug Fixes

<!-- none yet -->

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
