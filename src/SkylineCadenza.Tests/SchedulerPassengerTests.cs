using SkylineCadenza.Core.Scheduling;
using Xunit;

namespace SkylineCadenza.Tests;

/// <summary>
/// Passengers: precursors that already sit inside a scheduled slot's emitted
/// band and elute entirely inside its firing window, so the quadrupole
/// transmits them whether or not the assay asks for them.
///
/// The scheduler only ever evaluates a precursor as a DRIVER - something a slot
/// must grow to reach - and only against the slots that existed at its
/// protein's turn. The driver-side rules that shape a window (same charge per
/// slot, fragment clash, per-protein cap) therefore exclude precursors that are
/// co-isolated regardless, which removes nothing from the spectrum and only
/// stops the assay listing them. These tests pin the post-pass that recovers
/// them without touching the acquisition.
/// </summary>
public class SchedulerPassengerTests
{
    private static Candidate Make(
        string id, double mz, double rtStart, double rtStop, string group,
        double quantity, double[] fragments, int charge = 2)
    {
        var ions = new FragmentIon[fragments.Length];
        for (int i = 0; i < fragments.Length; i++)
            ions[i] = new FragmentIon(fragments[i], fragments.Length - i, 1);
        return new Candidate
        {
            PrecursorId = id,
            StrippedSequence = id,
            ModifiedSequence = id,
            PrecursorCharge = charge,
            PrecursorMz = mz,
            RtStart = rtStart,
            RtStop = rtStop,
            RtApex = (rtStart + rtStop) / 2.0,
            PrecursorQuantity = quantity,
            QValue = 0.001,
            ProteinQValue = 0.001,
            Proteotypic = 1,
            ProteinGroup = group,
            PeptideType = "unique",
            Fragments = ions,
            Top4Fragments = Candidate.DeriveTopMz(ions, 4),
            Run = "test",
        };
    }

    private static readonly double[] DriverFrags = { 100.0, 200.0, 300.0, 400.0 };
    private static readonly double[] CleanFrags = { 150.0, 250.0, 350.0, 450.0 };

    // "driver" opens a slot at 500.0 with band [499.65, 500.35] and firing
    // window [4.75, 5.65]. A rider at 500.2 is inside that band, so it can
    // never open a disjoint window of its own: it is either a member of the
    // driver's slot or nothing. Every case below turns on why it isn't a member.
    private static Candidate Driver() =>
        Make("driver", 500.0, 5.0, 5.4, "G1", 1e6, DriverFrags);

    private static ScheduleResult Run(params Candidate[] cands) =>
        Scheduler.Run(cands, new SchedulingParameters { CycleBudget = 100 });

    private static string[] Ids(ScheduleResult r, Candidate[] c, int[] idx) =>
        idx.Select(i => c[i].PrecursorId).ToArray();

    [Fact]
    public void Scheduler_Passenger_AdmitsPrecursorBlockedOnlyByChargeMismatch()
    {
        // Same band, co-eluting, no fragment clash: the ONLY thing stopping it
        // being a member is SameChargePerSlot. That rule exists so a slot has a
        // coherent collision energy, but a +3 precursor inside a +2 slot's band
        // is fragmented at that CE regardless.
        var cands = new[] { Driver(), Make("rider", 500.2, 5.0, 5.4, "G2", 1e5, CleanFrags, charge: 3) };

        var r = Run(cands);

        Assert.Equal(new[] { "driver" }, Ids(r, cands, r.ScheduledIndices));
        Assert.Equal(new[] { "rider" }, Ids(r, cands, r.PassengerIndices));
        Assert.False(r.PassengerClashesWithDriver[0]);
    }

    [Fact]
    public void Scheduler_Passenger_AdmitsClashingPrecursorAndFlagsIt()
    {
        // Fragment clash blocks membership but not co-isolation. Admitted, and
        // flagged so the ambiguous chromatogram is reviewable.
        var cands = new[] { Driver(), Make("rider", 500.2, 5.0, 5.4, "G2", 1e5, DriverFrags) };

        var r = Run(cands);

        Assert.Equal(new[] { "rider" }, Ids(r, cands, r.PassengerIndices));
        Assert.True(r.PassengerClashesWithDriver[0]);
    }

    [Fact]
    public void Scheduler_Passenger_IgnoresPerProteinCap()
    {
        // The cap limits how much budget one protein may spend. A passenger
        // spends none, so it must not be gated by it.
        var cands = new[] { Driver(), Make("rider", 500.2, 5.0, 5.4, "G1", 1e5, CleanFrags, charge: 3) };

        var r = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            MinPeptidesPerProtein = 1,
            MaxPeptidesPerProtein = 1,
        });

        Assert.Equal(new[] { "driver" }, Ids(r, cands, r.ScheduledIndices));
        Assert.Equal(new[] { "rider" }, Ids(r, cands, r.PassengerIndices));
    }

    [Fact]
    public void Scheduler_Passenger_RejectsPartialRtContainment()
    {
        // Peak runs past the end of the firing window, so the window is only on
        // for part of its elution. This is the strongest single predictor of a
        // ridealong not being re-identified, so it is a hard requirement.
        var cands = new[] { Driver(), Make("rider", 500.2, 5.5, 6.3, "G2", 1e5, CleanFrags, charge: 3) };

        var r = Run(cands);

        Assert.Equal(new[] { "driver" }, Ids(r, cands, r.ScheduledIndices));
        Assert.Empty(r.PassengerIndices);
    }

    [Fact]
    public void Scheduler_Passenger_RejectsPrecursorBelowShareFloor()
    {
        // 0.1% of the band's current: present, but with too little of the fill
        // to come back as a usable measurement.
        var cands = new[] { Driver(), Make("rider", 500.2, 5.0, 5.4, "G2", 1e3, CleanFrags, charge: 3) };

        var r = Run(cands);

        Assert.Empty(r.PassengerIndices);
    }

    [Fact]
    public void Scheduler_Passenger_DoesNotAlterSlotGeometry()
    {
        // The whole premise: a passenger is already inside the emitted band, so
        // admitting it must not move an edge. If it widened the slot it would
        // change the emitted isolation width and stop being free.
        var cands = new[] { Driver(), Make("rider", 500.2, 5.0, 5.4, "G2", 1e5, CleanFrags, charge: 3) };

        var r = Run(cands);

        var slot = Assert.Single(r.Slots);
        Assert.Equal(500.0, slot.MzMin, 9);
        Assert.Equal(500.0, slot.MzMax, 9);
        Assert.Equal(new[] { 0 }, slot.MemberIndices.ToArray());
        Assert.Equal(new[] { 1 }, slot.PassengerIndices.ToArray());
    }

    [Fact]
    public void Scheduler_Passenger_IsNeverAlsoADriver()
    {
        // Two precursors far enough apart to each drive their own slot, plus a
        // rider. Drivers and passengers must stay disjoint sets.
        var cands = new[]
        {
            Driver(),
            Make("other", 700.0, 5.0, 5.4, "G3", 8e5, new[] { 110.0, 210.0, 310.0, 410.0 }),
            Make("rider", 500.2, 5.0, 5.4, "G2", 1e5, CleanFrags, charge: 3),
        };

        var r = Run(cands);

        Assert.Equal(2, r.ScheduledIndices.Length);
        Assert.Empty(r.ScheduledIndices.Intersect(r.PassengerIndices));
    }

    [Fact]
    public void Scheduler_WithPassengers_FoldsRidersInAsScheduled()
    {
        // The analysis-side list. Same slots, same geometry - only the set of
        // precursors Skyline is told to extract changes.
        var cands = new[] { Driver(), Make("rider", 500.2, 5.0, 5.4, "G2", 1e5, CleanFrags, charge: 3) };

        var r = Run(cands);
        var combined = r.WithPassengers();

        Assert.Single(r.ScheduledIndices);
        Assert.Equal(2, combined.ScheduledIndices.Length);
        Assert.Equal(new[] { "driver", "rider" },
            Ids(combined, cands, combined.ScheduledIndices).OrderBy(x => x).ToArray());
        // Every rider is attributed to the slot it actually rides in.
        Assert.Equal(r.Slots[0].Id, combined.ScheduledSlotIds[^1]);
        Assert.Same(r.Slots, combined.Slots);
    }

    [Fact]
    public void Scheduler_PrmMode_AdmitsNoPassengers()
    {
        // PRM isolates one precursor per scan at a fixed quadrupole width. The
        // concept does not apply, and admitting riders would misrepresent what
        // the method targets.
        var cands = new[] { Driver(), Make("rider", 500.2, 5.0, 5.4, "G2", 1e5, CleanFrags, charge: 3) };

        var r = Scheduler.Run(cands, new SchedulingParameters
        {
            CycleBudget = 100,
            Mode = AcquisitionMode.Prm,
        });

        Assert.Empty(r.PassengerIndices);
    }
}
