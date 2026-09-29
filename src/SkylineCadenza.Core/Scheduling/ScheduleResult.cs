namespace SkylineCadenza.Core.Scheduling;

/// <summary>
/// Output of <see cref="Scheduler.Run(IReadOnlyList{Candidate}, SchedulingParameters, System.Threading.CancellationToken)"/>.
/// </summary>
public sealed class ScheduleResult
{
    /// <summary>Indexes into the input <c>candidates</c> list of every scheduled precursor.</summary>
    public required int[] ScheduledIndices { get; init; }

    /// <summary>Slot id assigned to each scheduled precursor (parallel array to <see cref="ScheduledIndices"/>).</summary>
    public required int[] ScheduledSlotIds { get; init; }

    /// <summary>All slots created during scheduling.</summary>
    public required Slot[] Slots { get; init; }

    /// <summary>Centers of the RT-occupancy grid (min).</summary>
    public required double[] RtGrid { get; init; }

    /// <summary>Number of concurrent slots active in each RT bin.</summary>
    public required int[] SlotCountCurve { get; init; }

    /// <summary>Distinct protein groups represented in the schedule.</summary>
    public required int ProteinGroupsCovered { get; init; }

    /// <summary>
    /// Indexes of precursors admitted as passengers (see
    /// <see cref="Slot.PassengerIndices"/>). Disjoint from
    /// <see cref="ScheduledIndices"/>: a precursor is either a driver or a
    /// passenger, never both.
    /// </summary>
    public int[] PassengerIndices { get; init; } = Array.Empty<int>();

    /// <summary>Slot id each passenger rides in (parallel to <see cref="PassengerIndices"/>).</summary>
    public int[] PassengerSlotIds { get; init; } = Array.Empty<int>();

    /// <summary>
    /// Whether each passenger shares a fragment m/z (within
    /// <see cref="SchedulingParameters.FragmentTolDa"/>) with one of its
    /// slot's drivers. Not a rejection reason - measured against the acquired
    /// 08/03 MTM run, clashing passengers were re-identified just as often as
    /// clean ones - but it flags the ones whose extracted chromatograms may be
    /// ambiguous, so they can be reviewed.
    /// </summary>
    public bool[] PassengerClashesWithDriver { get; init; } = Array.Empty<bool>();

    /// <summary>
    /// One outcome per <see cref="SchedulingParameters.ForcedPeptides"/> entry,
    /// in list order. Empty when no forced list was given.
    /// </summary>
    public ForcedPeptideOutcome[] ForcedOutcomes { get; init; } = Array.Empty<ForcedPeptideOutcome>();

    /// <summary>
    /// The same schedule with passengers folded in as scheduled precursors, for
    /// the analysis-side target list. The acquisition is identical either way:
    /// passengers are already inside the emitted bands, so no window, width,
    /// collision energy, or cycle-budget value changes. Use the unmodified
    /// result for the instrument method and this one for the Skyline document.
    /// </summary>
    public ScheduleResult WithPassengers()
    {
        if (PassengerIndices.Length == 0) return this;
        var idx = new int[ScheduledIndices.Length + PassengerIndices.Length];
        var sid = new int[idx.Length];
        ScheduledIndices.CopyTo(idx, 0);
        ScheduledSlotIds.CopyTo(sid, 0);
        PassengerIndices.CopyTo(idx, ScheduledIndices.Length);
        PassengerSlotIds.CopyTo(sid, ScheduledSlotIds.Length);
        return new ScheduleResult
        {
            ScheduledIndices = idx,
            ScheduledSlotIds = sid,
            Slots = Slots,
            RtGrid = RtGrid,
            SlotCountCurve = SlotCountCurve,
            ProteinGroupsCovered = ProteinGroupsCovered,
            PassengerIndices = PassengerIndices,
            PassengerSlotIds = PassengerSlotIds,
            PassengerClashesWithDriver = PassengerClashesWithDriver,
            ForcedOutcomes = ForcedOutcomes,
        };
    }
}
