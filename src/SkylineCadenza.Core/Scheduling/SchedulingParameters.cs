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

    /// <summary>
    /// Rank cap: only the top N candidates of each protein group are ELIGIBLE
    /// for scheduling. 0 (the default) means no cap, which is what keeps every
    /// existing schedule and the golden fixtures unchanged.
    ///
    /// The cap counts PRECURSORS, not distinct sequences: a peptide present at
    /// two charge states consumes two of the N. Ranking is whatever
    /// <see cref="PeptideRanking"/> selects (quantity, q-value, or predicted
    /// score), applied after the unique-before-razor tier, so the knob works
    /// for every peptide-priority option. It is applied PER SOURCE TIER over
    /// that tier's own ranking metric, so in DIA-NN + Carafe backfill mode the
    /// confirmed tier keeps its top N by evidence and the predicted tier keeps
    /// its top N by predicted score, with neither tier's cap perturbing the
    /// other.
    ///
    /// This truncates the candidate pool, so it also feeds the protein-ordering
    /// metrics (summed intensity / best q-value are computed over the capped
    /// queue) and the smallest-group tiebreaker. Distinct from
    /// <see cref="MaxPeptidesPerProtein"/>, which caps how many of the eligible
    /// peptides actually get scheduled.
    /// </summary>
    public int PeptideRankCap { get; init; } = 0;

    /// <summary>
    /// Minimum share of the co-isolated ion current that EVERY measured target
    /// in a slot must hold while it elutes, as a fraction in [0, 1). A
    /// placement is refused when it would leave any member below this share,
    /// so a target is never scheduled into a window where competing ions would
    /// take so much of the fill that it cannot be quantified.
    ///
    /// The denominator is every pooled precursor whose m/z falls in the slot's
    /// emitted band (<c>MzMin - PrmIsolationWidthTh / 2</c> to
    /// <c>MzMax + PrmIsolationWidthTh / 2</c>), whether or not it is itself
    /// scheduled, each weighted by the fraction of ITS peak that co-elutes
    /// with the member being judged. Both scopes matter: on a real HeLa MTM
    /// acquisition roughly two thirds of the ion current inside the target
    /// windows came from precursors that were never targets, so a
    /// members-only check measures almost none of the competition; and
    /// counting a neighbour's whole area on any overlap at all overstated the
    /// competition by more than 3x for 57% of targets, scoring a huge
    /// neighbour that clips the edge of a peak the same as one sitting on it.
    ///
    /// 0 (the default) disables the check, which keeps existing schedules and
    /// the golden fixtures unchanged; the app defaults it to 0.02 (2%).
    /// Precursors with no measured intensity (quantity &lt;= 0, i.e. Carafe
    /// predictions) are invisible to it: they neither count as background nor
    /// get protected as targets.
    ///
    /// The estimate counts only IDENTIFIED precursors, so it is a lower bound
    /// on the true co-isolated current; unidentified background sits on top of
    /// it and the threshold has to absorb that.
    /// </summary>
    public double MinTargetIonShare { get; init; } = 0.0;

    /// <summary>
    /// Fragment-level co-isolation gate (MTM only). When on, a placement is
    /// refused if any member's top-4 fragment falls within
    /// <see cref="FragmentTolDa"/> of a top-4 fragment of ANY other identified
    /// precursor that co-elutes with it inside the slot's emitted band
    /// (<c>MzMin - PrmIsolationWidthTh / 2</c> to
    /// <c>MzMax + PrmIsolationWidthTh / 2</c>). This covers other drivers AND
    /// never-scheduled background, which the driver-vs-driver clash check
    /// cannot see. Every top-4 fragment of a target must be clean, where a
    /// fragment counts as interfered once the overlap-weighted co-isolated
    /// signal on it reaches 1% of its own signal (a fixed tolerance, not a
    /// knob), so a far weaker background precursor does not block a target.
    /// Only precursors with a measured quantity are counted as interferers.
    ///
    /// False (the default) keeps existing schedules and the golden fixtures
    /// unchanged; the app always turns it on.
    /// </summary>
    public bool BackgroundFragmentClash { get; init; } = false;

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

    /// <summary>Convenience accessor: <see cref="FiringPadSec"/> in minutes.</summary>
    public double FiringPadMin => FiringPadSec / 60.0;
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
