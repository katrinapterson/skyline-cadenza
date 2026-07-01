# Skyline Cadenza Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on `main`. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and `Version` in `src/SkylineCadenza.App/tool-inf/info.properties` is updated to match.

## New Features

- Added a separate **Carafe detectability predictions** input (multi-file / folder) plus a **detectability cutoff** slider (default 0.5) for the Carafe path. This is independent of the predicted peptide-ranking overlay, so a user can rank with Bromo and filter by Carafe detectability simultaneously. Precursors below the cutoff are dropped before scheduling (peptides with no prediction are kept); the slider rebuilds the candidate pool and reschedules live. Empirical calibration: ~0.5 removes peptides ~90% unlikely to be observed while keeping most real detections.

- Added a configurable **Predicted peak half-width** slider (seconds) for the Carafe path, replacing the hardcoded 0.10 min synthesised peak. Adjusting it rebuilds the Carafe candidate pool and reschedules live; ignored on the DIA-NN / Skyline paths, which use measured peak boundaries.
- Added a predicted peptide-ranking overlay for Carafe-only mode. Supply a Bromo `peptide_rankings.tsv` or Carafe `*_detectability_pred.tsv` (multi-select several shard files, a folder of shards, or one path per line) via the new "Predicted peptide ranking" input. The score is matched onto the Carafe candidate pool by (peptide, charge) and drives a real within-protein peptide priority (`PeptidePriority.PredictedRanking`), which auto-selects when an overlay is loaded. Replaces the no-signal fragment-intensity-sum stand-in documented in `docs/ingest-sources.md`. Within-protein only; it does not affect cross-protein protein ranking.

## Bug Fixes

<!-- none yet -->

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
