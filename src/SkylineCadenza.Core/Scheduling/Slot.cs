namespace SkylineCadenza.Core.Scheduling;

/// <summary>
/// An MTM acquisition slot. Holds one or more co-eluting precursors that
/// share a single isolation window definition.
/// </summary>
/// <remarks>
/// <para>
/// <see cref="RtStart"/> / <see cref="RtStop"/> are the padded firing window
/// (union of member padded ranges). This is what the budget counter tracks -
/// the time the window is on.
/// </para>
/// <para>
/// <see cref="CoStart"/> / <see cref="CoStop"/> are the intersection of all
/// members' UNPADDED peak ranges. This is when every member is actually
/// co-eluting; multiplexed firing produces useful deconvolvable data only
/// during this interval. A non-empty intersection is required for any join.
/// </para>
/// </remarks>
public sealed class Slot
{
    public int Id { get; set; }

    public List<int> MemberIndices { get; } = new();

    public double MzMin { get; set; }
    public double MzMax { get; set; }

    /// <summary>
    /// Hard floor on <see cref="MzMin"/>. Set from the nearest co-eluting,
    /// clashing precursor BELOW this slot in the candidate pool: the emitted
    /// band reaches <c>PrmIsolationWidthTh / 2</c> past the edge member, so the
    /// bound sits half a quadrupole width inside the barrier. Unbounded in
    /// PRM mode, where a slot holds one precursor and never grows.
    /// </summary>
    public double GrowLo { get; set; } = double.NegativeInfinity;

    /// <summary>Hard ceiling on <see cref="MzMax"/>. See <see cref="GrowLo"/>.</summary>
    public double GrowHi { get; set; } = double.PositiveInfinity;

    /// <summary>Padded firing window start (min).</summary>
    public double RtStart { get; set; }

    /// <summary>Padded firing window stop (min).</summary>
    public double RtStop { get; set; }

    /// <summary>Unpadded co-elution intersection start (min).</summary>
    public double CoStart { get; set; }

    /// <summary>Unpadded co-elution intersection stop (min).</summary>
    public double CoStop { get; set; }

    /// <summary>
    /// Precursors that ride along in this slot: they fall inside the emitted
    /// band and elute entirely within the firing window, so the quadrupole
    /// already transmits them and they are fragmented in this scan whether or
    /// not the assay asks for them.
    /// </summary>
    /// <remarks>
    /// Passengers are NOT members. <see cref="MemberIndices"/> are the drivers:
    /// they set <see cref="MzMin"/> / <see cref="MzMax"/>, the emitted
    /// isolation width, the collision energy, and the cycle budget. Passengers
    /// are admitted in a post-pass once the schedule is frozen and change none
    /// of those, which is why they cost nothing and why the driver-side rules
    /// that exist to shape a window (same charge per slot, fragment clash,
    /// per-protein cap) do not apply to them.
    /// </remarks>
    public List<int> PassengerIndices { get; } = new();

    /// <summary>
    /// Union of all member top-4 fragments, sorted ascending. New joins must
    /// not clash with any value in here within
    /// <see cref="SchedulingParameters.FragmentTolDa"/>.
    /// </summary>
    public double[] Fragments { get; set; } = Array.Empty<double>();
}
