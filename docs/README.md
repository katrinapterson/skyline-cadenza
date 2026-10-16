# Skyline Cadenza Documentation

User and developer documentation for Skyline Cadenza. For repo orientation
(what the project is, how it's organised, how to build and test it), see
`CLAUDE.md` at the repo root; this folder covers the algorithmic and
integration details that don't belong in a development-context file.

## Contents

- [`skyline-integration.md`](skyline-integration.md): how Cadenza connects
  to a running Skyline document, what it reads from the document, and what
  it writes back when the user clicks "Update Skyline document". Covers
  the JSON-RPC pipe, the BLIB writer, the `RunCommand`-driven settings
  configurator, and the transition-list import path.

- [`coverage-objectives.md`](coverage-objectives.md): the three greedy
  scheduling objectives (Balanced, Maximize Proteins, Maximize Peptides),
  what each one does at the cover-pass and load-up stages, and when to
  pick which.

- [`settings.md`](settings.md): every user-tunable knob exposed by the UI
  and the underlying `SchedulingParameters` record, with sensible default
  values, units, and notes on which settings interact.

- [`ingest-sources.md`](ingest-sources.md): the three peptide-pool input
  paths (DIA-NN report, Carafe predicted library, running Skyline
  document), what they contribute to a candidate, and the answer to "how
  does a Carafe-only run get peptide rankings when Carafe has no
  abundance information?".

- [`scheduler-window-scoping.md`](scheduler-window-scoping.md): known bug -
  `TrySchedule`'s fragment-clash check is scoped to a slot's members, but
  the isolation window physically co-isolates everything in its m/z band, so
  expanding a window can silently re-include a precursor that was rejected
  for clashing (the "illusory exclusion"). Motivates the window-scoped
  interference redesign. Guarded by a skipped regression test.

- [`automation.md`](automation.md): driving Cadenza from another app
  (AutoTargeting). Covers the headless `SkylineCadenzaCli export-thermo`
  verb, its PARAMS.json fields and app defaults, and the GUI's standalone
  `--report` / `--export-dir` launch options.

## Screenshots

`screenshots/` carries the figures referenced from the release notes and
the docs above (UI, MTM slot layout, protein-coverage plot).
