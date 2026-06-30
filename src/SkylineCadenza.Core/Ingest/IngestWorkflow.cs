namespace SkylineCadenza.Core.Ingest;

/// <summary>
/// Which file-based candidate-pool workflow the "Load + schedule" command runs.
/// Selected explicitly in the UI rather than inferred from which paths are set.
/// (The "Load from Skyline document" path is a separate entry point and is not
/// governed by this switch.)
/// </summary>
public enum IngestWorkflow
{
    /// <summary>
    /// DIA-NN report drives the candidate pool: measured RT boundaries,
    /// intensity, and q-values. A Carafe library is optional on this path and,
    /// when supplied, only overrides per-row fragments.
    /// </summary>
    Diann,

    /// <summary>
    /// Carafe predicted spectral library + a target protein list drive the
    /// pool. There is no measured peak; the firing window is synthesised from
    /// the predicted apex by the user-set peak half-width. An optional
    /// predicted-ranking overlay sets within-protein peptide priority.
    /// </summary>
    Carafe,
}
