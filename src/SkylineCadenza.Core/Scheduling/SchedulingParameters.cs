namespace SkylineCadenza.Core.Scheduling;

/// <summary>
/// User-tunable scheduler knobs. Defaults match the targeted-modeling
/// notebook's "MTM at budget 100" configuration.
/// </summary>
public sealed record SchedulingParameters
{
    /// <summary>PRM (1 precursor per slot) or MTM (multiplex).</summary>
    public AcquisitionMode Mode { get; init; } = AcquisitionMode.Mtm;

    /// <summary>Maximum MTM isolation window width (Th). Ignored when <see cref="Mode"/> is PRM.</summary>
    public double IsolationWindowTh { get; init; } = 3.0;

    /// <summary>
    /// Quadrupole isolation window width (Th) the instrument fires around
    /// each precursor m/z. Used in three places:
    /// <list type="bullet">
    /// <item>PRM mode: the isolation width for every scheduled precursor.</item>
    /// <item>MTM solo slots: the instrument isolation written to the
    /// Thermo CSV for slots holding a single precursor.</item>
    /// <item>MTM slot-edge rule: every member's quadrupole window must
    /// fit fully inside the slot's nominal <see cref="IsolationWindowTh"/>,
    /// so the scheduler keeps every precursor center at least
    /// (PrmIsolationWidthTh / 2) from each slot edge.</item>
    /// </list>
    /// </summary>
    public double PrmIsolationWidthTh { get; init; } = 0.7;

    /// <summary>Fragment m/z tolerance for the non-clash check (Da).</summary>
    public double FragmentTolDa { get; init; } = 0.5;

    /// <summary>Max concurrent acquisition slots per cycle.</summary>
    public int CycleBudget { get; init; } = 100;

    /// <summary>
    /// Seconds of padding added to each side of each peak's retention range,
    /// added on top of the measured peak boundaries. For DIA-NN / Skyline-BLIB
    /// ingest, the boundaries already include per-replicate peak-shape variance
    /// (Cadenza takes the union of every replicate's <c>startTime</c> /
    /// <c>endTime</c>), so this padding handles unmodelled drift like
    /// gradient drift between the source acquisitions and the new run, not
    /// the peak shape itself.
    /// </summary>
    public double FiringPadSec { get; init; } = 15.0;

    /// <summary>Width of RT-occupancy bins (min). 0.05 = 3 s.</summary>
    public double RtBinMin { get; init; } = 0.05;

    /// <summary>DIA-NN <c>Q.Value</c> cutoff used at the candidate-build stage.</summary>
    public double QValueCutoff { get; init; } = 0.01;

    /// <summary>
    /// User-supplied set of protein accessions to focus the assay on. Empty
    /// = no filter (the scheduler considers every protein the report sees).
    /// </summary>
    public IReadOnlySet<string> TargetProteins { get; init; } = new HashSet<string>();

    /// <summary>
    /// What to do with proteins that AREN'T in <see cref="TargetProteins"/>.
    /// Has no effect if the target set is empty.
    /// </summary>
    public TargetListMode TargetMode { get; init; } = TargetListMode.FirstThenFill;

    /// <summary>Minimum peptides scheduled per protein group (best effort).</summary>
    public int MinPeptidesPerProtein { get; init; } = 1;

    /// <summary>Maximum peptides scheduled per protein group.</summary>
    public int MaxPeptidesPerProtein { get; init; } = 5;

    /// <summary>How to order the protein groups going into the assay.</summary>
    public ProteinPriority ProteinRanking { get; init; } = ProteinPriority.SummedIntensity;

    /// <summary>How to order peptides within each protein group.</summary>
    public PeptidePriority PeptideRanking { get; init; } = PeptidePriority.PrecursorQuantity;

    /// <summary>
    /// How precursor charge is handled when joining MTM slots. Ignored for
    /// PRM (always 1 precursor / slot).
    /// </summary>
    public ChargeHandling ChargeHandling { get; init; } = ChargeHandling.SameChargePerSlot;

    /// <summary>
    /// Normalized collision energy emitted in the Thermo inclusion CSV.
    /// Not used by the scheduler; only by <c>ThermoCsvWriter</c>.
    /// </summary>
    public double NormalizedCollisionEnergy { get; init; } = 28.0;

    /// <summary>
    /// What the scheduler optimises for after meeting the
    /// load-balanced RT-budget constraint. See
    /// <see cref="CoverageObjective"/> for the precise behaviour of
    /// each mode. Determines both the cover-pass strategy
    /// (reactive vs look-ahead) and the load-up cap, so no separate
    /// <c>EnableLoadBalancing</c> / <c>RtAwareCoverSelection</c>
    /// toggles are needed.
    /// </summary>
    public CoverageObjective Objective { get; init; } = CoverageObjective.MaximizeProteins;

    /// <summary>
    /// How MTM isolation windows are allowed to grow relative to each other
    /// and to precursors they would sweep over. See <see cref="WindowPacking"/>.
    /// Defaults to <see cref="WindowPacking.Legacy"/> so existing schedules and
    /// the golden fixtures are unchanged.
    /// </summary>
    public WindowPacking Packing { get; init; } = WindowPacking.Legacy;

    /// <summary>Convenience accessor: <see cref="FiringPadSec"/> in minutes.</summary>
    public double FiringPadMin => FiringPadSec / 60.0;
}

/// <summary>
/// Controls whether a slot's isolation window may expand across precursors it
/// would co-isolate, and whether two concurrent windows may overlap in m/z.
/// </summary>
public enum WindowPacking
{
    /// <summary>
    /// Historical behaviour: each slot grows to cover its own members with no
    /// awareness of its neighbours. A merge that widens the band can sweep it
    /// over a precursor that was rejected for clashing (the "illusory
    /// exclusion" in <c>docs/scheduler-window-scoping.md</c>), and two
    /// concurrent windows may overlap, double-sampling the precursors between
    /// them.
    /// </summary>
    Legacy,

    /// <summary>
    /// Slot-aware packing. Every slot carries growth bounds
    /// (<see cref="Slot.GrowLo"/> / <see cref="Slot.GrowHi"/>) derived from the
    /// nearest co-eluting, clashing precursor in the candidate pool on each
    /// side, so a window can never expand across an interferer. Expansions and
    /// new slots are additionally rejected when the resulting band would
    /// overlap a concurrently-firing slot's band, so windows never grow into
    /// one another. Protein-priority order and the per-protein fallback are
    /// unchanged; the cost is paid in cycle budget, not coverage logic.
    /// </summary>
    Disjoint,
}

public enum AcquisitionMode
{
    Prm,
    Mtm,
}

public enum TargetListMode
{
    /// <summary>Only proteins in <c>TargetProteins</c> are eligible.</summary>
    Exclusive,
    /// <summary>Target proteins first, then fill with off-list proteins.</summary>
    FirstThenFill,
}

public enum ProteinPriority
{
    /// <summary>Sum of peptide-precursor quantities per group (DIA-NN intensity).</summary>
    SummedIntensity,
    /// <summary>Best (smallest) per-group q-value (<c>PG.Q.Value</c>).</summary>
    ProteinQValue,
    /// <summary>Order proteins by the user's input list (insertion order).</summary>
    ProvidedListOrder,
}

public enum PeptidePriority
{
    /// <summary>DIA-NN <c>Precursor.Quantity</c> (real measured intensity).</summary>
    PrecursorQuantity,
    /// <summary>DIA-NN <c>Q.Value</c> (peptide-precursor q-value, ascending).</summary>
    QValue,
    /// <summary>
    /// Predicted within-protein ranking (<see cref="Candidate.PredictedPriority"/>,
    /// descending; candidates with no predicted score sort last). Sourced from a
    /// Bromo <c>mean_pred_score</c> or Carafe <c>detectability_pred</c> overlay in
    /// Carafe-only mode, where it auto-selects as the active peptide priority when
    /// an overlay file is loaded.
    /// </summary>
    PredictedRanking,
}

public enum CoverageObjective
{
    /// <summary>
    /// Exact published Stellar-webinar algorithm. Reactive cover pass:
    /// walk each protein's score-sorted queue and take the first peptide
    /// that fits the RT budget; on a per-bin budget hit, fall back to
    /// the next peptide in score order (webinar slide 9). Round-robin
    /// load-up to <see cref="SchedulingParameters.MaxPeptidesPerProtein"/>.
    /// Use this when you want a reproducible baseline matching the
    /// Remes / MacCoss webinar; otherwise prefer
    /// <see cref="MaximizeProteins"/>.
    /// </summary>
    Balanced,

    /// <summary>
    /// Maximise the number of protein groups with at least one peptide
    /// scheduled. Look-ahead cover pass: prefers peptides that join an
    /// existing slot (no new slot opened, no cycle budget consumed),
    /// falling back to the least-saturated peptide when no joinable
    /// option exists. After every coverable protein has its first
    /// peptide, the round-robin load-up fills up to
    /// <see cref="SchedulingParameters.MaxPeptidesPerProtein"/> per
    /// group using whatever cycle budget remains. This is Cadenza's
    /// default.
    /// </summary>
    MaximizeProteins,

    /// <summary>
    /// Maximise the total number of scheduled peptides. Look-ahead
    /// cover pass (same as <see cref="MaximizeProteins"/> minus the
    /// joinable-first preference). Load-up uncapped (the per-group max
    /// is ignored): the round-robin continues until either the cycle
    /// budget saturates or
    /// every group has exhausted its candidate queue. Round-robin order
    /// is preserved so a protein never receives its k+1th peptide before
    /// every other covered protein has had a chance at its kth.
    /// </summary>
    MaximizePeptides,
}

public enum ChargeHandling
{
    /// <summary>
    /// MTM slots may only contain precursors of identical
    /// <see cref="Candidate.PrecursorCharge"/>. Cleaner CE assignment and
    /// simpler downstream interpretation; mild coverage cost.
    /// </summary>
    SameChargePerSlot,
    /// <summary>
    /// MTM slots may contain mixed charges. The Thermo CSV writer reports
    /// each slot's majority charge (ties go to the lower z, since +2 is the
    /// canonical tryptic default).
    /// </summary>
    AllowMixed,
}
