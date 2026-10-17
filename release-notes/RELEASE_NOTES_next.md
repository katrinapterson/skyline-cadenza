# Skyline Cadenza Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on `main`. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and `Version` in `src/SkylineCadenza.App/tool-inf/info.properties` is updated to match.

## New Features

- **Reorganized the parameter panel around an explicit data source.** A new "Data source" choice (Empirical, Empirical + predicted backfill, Predicted only, Skyline document) replaces guessing the mode from which file boxes are filled and the "Backfill unused schedule time with Carafe predictions" checkbox. Each source shows only the inputs it reads: the Bromo ranking, detectability, and predicted peak half-width controls appear only when predictions are loaded, and the Q value cutoff only when there is empirical evidence. "Load + schedule" now also loads from the Skyline document. The remaining settings are grouped into collapsible Targets, Acquisition, Selection strategy, Cycle budget, and Output sections. MTM-only controls are hidden in PRM mode, Normalized CE moved to Acquisition, and the PredictedRanking peptide ranking is offered only in predicted-only mode, where it has an effect.
- **Cadenza can be driven by another app.** `SkylineCadenzaCli export-thermo REPORT.parquet PARAMS.json|- OUTPUT.csv` writes the same Thermo scheduled inclusion CSV as **Export Thermo CSV...**, using the app's default settings unless PARAMS.json sets them. `SkylineCadenza.exe --report <report.parquet> [--export-dir <folder>]` opens Cadenza standalone without Skyline, loads and schedules the report, and opens **Export Thermo CSV...** in the given folder. Launching from Skyline is unchanged. See docs/automation.md.
- **Another app can see what a design scheduled and choose targets by Cadenza's own names.** `export-thermo` takes an optional fourth argument, SUMMARY.json, describing what was scheduled: each target-list entry's matched protein groups and scheduled precursors, each forced peptide's outcome, and the peptides per protein. `SkylineCadenzaCli describe-diann REPORT.parquet PARAMS.json|- OUTPUT.json` lists the protein groups (with gene names) and precursors the scheduler would choose from. See docs/automation.md.
- **The command line can write the assay as a Skyline document.** `export-thermo ... --skyline OUT.sky` writes a document with the assay library, the settings **Push to Skyline** applies and the targets, through SkylineCmd, so the targeted runs can be imported without opening Skyline. See docs/automation.md.

## Bug Fixes

## Performance

## Breaking Changes
