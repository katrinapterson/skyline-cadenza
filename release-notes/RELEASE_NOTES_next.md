# Skyline Cadenza Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on `main`. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and `Version` in `src/SkylineCadenza.App/tool-inf/info.properties` is updated to match.

## New Features

- Added a predicted peptide-ranking overlay for Carafe-only mode. Supply a Bromo `peptide_rankings.tsv` or Carafe `*_detectability_pred.tsv` (multi-select several shard files, a folder of shards, or one path per line) via the new "Predicted peptide ranking" input. The score is matched onto the Carafe candidate pool by (peptide, charge) and drives a real within-protein peptide priority (`PeptidePriority.PredictedRanking`), which auto-selects when an overlay is loaded. Replaces the no-signal fragment-intensity-sum stand-in documented in `docs/ingest-sources.md`. Within-protein only; it does not affect cross-protein protein ranking.

## Bug Fixes

<!-- none yet -->

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
