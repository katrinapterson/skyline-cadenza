#nullable enable
using SkylineCadenza.Core.Scheduling;
using Xunit;

namespace SkylineCadenza.Tests;

/// <summary>
/// Regression guard: the schedule must depend only on the CONTENT of the
/// candidate pool, never on the order the candidates arrive in.
/// </summary>
/// <remarks>
/// This is not hypothetical. <c>CarafeLibraryReader</c> materialises its
/// candidate list by enumerating a <c>Dictionary&lt;FragmentKey, ...&gt;</c>,
/// and .NET randomizes string hashing per process - so the ingest order
/// genuinely changes on every app launch. Before the peptide comparator was
/// made a total order, any pool with tied ranking scores produced a DIFFERENT
/// schedule from identical inputs across runs (observed jaccard 0.58-0.74),
/// which silently broke reproducibility between an exported target list and
/// the assay document built from it.
/// </remarks>
public class SchedulerDeterminismTests
{
    private static Candidate Make(string id, double mz, double rt, string group, double qty)
    {
        var frags = new[] { 100.0, 200.0, 300.0, 400.0 };
        var ions = new FragmentIon[frags.Length];
        for (int i = 0; i < frags.Length; i++) ions[i] = new FragmentIon(frags[i], frags.Length - i, 1);
        return new Candidate
        {
            PrecursorId = id, StrippedSequence = id, ModifiedSequence = id,
            PrecursorCharge = 2, PrecursorMz = mz,
            RtStart = rt, RtStop = rt + 0.4, RtApex = rt + 0.2,
            PrecursorQuantity = qty, QValue = 0.001, ProteinQValue = 0.001,
            Proteotypic = 1, ProteinGroup = group, PeptideType = "unique",
            Fragments = ions, Top4Fragments = Candidate.DeriveTopMz(ions, 4), Run = "t",
        };
    }

    /// <summary>
    /// Pool with heavy co-elution / m-z competition (so the cycle budget
    /// actually binds) and only three distinct quantities, so most peptides
    /// TIE on the ranking metric - the condition that used to expose the bug.
    /// </summary>
    private static List<Candidate> TiedScorePool()
    {
        var c = new List<Candidate>();
        for (int p = 0; p < 150; p++)
            for (int j = 0; j < 6; j++)
                c.Add(Make($"P{p}_pep{j}", 500.0 + (p % 12) * 0.4, 5.0 + (p % 30) * 0.5,
                           $"G{p}", 1e6 * (j % 3 + 1)));
        return c;
    }

    private static List<Candidate> Shuffled(List<Candidate> src, int seed)
    {
        var s = new List<Candidate>(src);
        var rnd = new Random(seed);
        for (int i = s.Count - 1; i > 0; i--) { int j = rnd.Next(i + 1); (s[i], s[j]) = (s[j], s[i]); }
        return s;
    }

    private static HashSet<string> Schedule(List<Candidate> pool, CoverageObjective objective,
                                            WindowPacking packing = WindowPacking.Legacy)
    {
        var r = Scheduler.Run(pool, new SchedulingParameters
        {
            CycleBudget = 8, Objective = objective, Packing = packing,
        });
        return r.ScheduledIndices.Select(i => pool[i].PrecursorId).ToHashSet();
    }

    [Theory]
    [InlineData(CoverageObjective.Balanced)]
    [InlineData(CoverageObjective.MaximizeProteins)]
    [InlineData(CoverageObjective.MaximizePeptides)]
    public void Scheduler_TiedScores_ScheduleIsIndependentOfCandidateOrder(CoverageObjective objective)
    {
        var pool = TiedScorePool();
        var expected = Schedule(pool, objective);

        foreach (int seed in new[] { 1, 42, 7919 })
        {
            var actual = Schedule(Shuffled(pool, seed), objective);
            Assert.True(expected.SetEquals(actual),
                $"{objective}: shuffling the candidate list (seed {seed}) changed the schedule " +
                $"({expected.Intersect(actual).Count()} of {expected.Count} shared). The peptide " +
                $"comparator must be a total order.");
        }
    }

    /// <summary>
    /// Same guarantee under <see cref="WindowPacking.Disjoint"/>. Slot-aware
    /// packing adds two order-sensitive-looking mechanisms - growth bounds and
    /// the neighbour-band probe - so it gets its own coverage. Bounds are
    /// sourced from the static candidate POOL (not from the slots built so
    /// far), and only ever tighten, which is what keeps them independent of
    /// the order proteins happen to be processed in.
    /// </summary>
    [Theory]
    [InlineData(CoverageObjective.Balanced)]
    [InlineData(CoverageObjective.MaximizeProteins)]
    [InlineData(CoverageObjective.MaximizePeptides)]
    public void Scheduler_DisjointPacking_ScheduleIsIndependentOfCandidateOrder(CoverageObjective objective)
    {
        var pool = TiedScorePool();
        var expected = Schedule(pool, objective, WindowPacking.Disjoint);
        Assert.NotEmpty(expected);

        foreach (int seed in new[] { 1, 42, 7919 })
        {
            var actual = Schedule(Shuffled(pool, seed), objective, WindowPacking.Disjoint);
            Assert.True(expected.SetEquals(actual),
                $"{objective} + Disjoint: shuffling the candidate list (seed {seed}) changed the " +
                $"schedule ({expected.Intersect(actual).Count()} of {expected.Count} shared).");
        }
    }
}
