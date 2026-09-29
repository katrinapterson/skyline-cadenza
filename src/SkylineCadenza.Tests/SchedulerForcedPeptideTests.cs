using SkylineCadenza.Core.Scheduling;
using Xunit;

namespace SkylineCadenza.Tests;

public class SchedulerForcedPeptideTests
{
    private static Candidate Make(
        string seq, int charge, double mz, double rtStart, double rtStop, string group,
        double quantity = 1e6, double[]? fragments = null, string? modSeq = null)
    {
        var mzList = fragments ?? new[] { 100.0, 200.0, 300.0, 400.0 };
        var ions = new FragmentIon[mzList.Length];
        for (int i = 0; i < mzList.Length; i++)
            ions[i] = new FragmentIon(mzList[i], mzList.Length - i, 1);
        return new Candidate
        {
            PrecursorId = (modSeq ?? seq) + charge,
            StrippedSequence = seq,
            ModifiedSequence = modSeq ?? seq,
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

    private static SchedulingParameters Params(string forced, int budget = 100, double share = 0.0) => new()
    {
        CycleBudget = budget,
        MinTargetIonShare = share,
        BackgroundFragmentClash = true,
        ForcedPeptides = ForcedPeptides.Parse(forced, out _),
    };

    private static HashSet<string> Scheduled(ScheduleResult r, IReadOnlyList<Candidate> c) =>
        r.ScheduledIndices.Select(i => c[i].PrecursorId).ToHashSet();

    [Fact]
    public void Parse_AcceptsChargeForms_AndSkipsHeadersAndComments()
    {
        var specs = ForcedPeptides.Parse(
            "Peptide,Charge\n# validation set\nPEPTIDEK\nPEPTIDER+2\nPEPTIDES/3\nPEPTIDET++\nPEPTIDEV,2\nPEPTIDEK\n",
            out var invalid);

        Assert.Equal(new[] { "Peptide,Charge" }, invalid);
        Assert.Equal(5, specs.Count);
        Assert.Null(specs[0].Charge);
        Assert.Equal(new int?[] { 2, 3, 2, 2 }, specs.Skip(1).Select(s => s.Charge));
    }

    [Fact]
    public void Matches_ModifiedEntryAcrossSyntaxes_UnmodifiedEntryMatchesAnyForm()
    {
        var c = Make("PEPCTIDEK", 2, 500.0, 5.0, 5.4, "G1", modSeq: "PEPC(UniMod:4)TIDEK");

        Assert.True(ForcedPeptides.Matches(ForcedPeptides.Parse("PEPC[+57.02]TIDEK", out _)[0], c));
        Assert.True(ForcedPeptides.Matches(ForcedPeptides.Parse("_PEPC[UniMod:4]TIDEK_", out _)[0], c));
        Assert.True(ForcedPeptides.Matches(ForcedPeptides.Parse("PEPCTIDEK", out _)[0], c));
        Assert.False(ForcedPeptides.Matches(ForcedPeptides.Parse("PEPC[+16.00]TIDEK", out _)[0], c));
        Assert.False(ForcedPeptides.Matches(ForcedPeptides.Parse("PEPCTIDEK+3", out _)[0], c));
    }

    [Fact]
    public void Scheduler_ForcedPeptide_IsScheduledEvenWhenOutrankedAndOverPerProteinMax()
    {
        // One protein, max one peptide. Ranking would pick TOPK; forcing LOWK
        // puts LOWK in, and the fill cannot add TOPK past the max.
        var cands = new[]
        {
            Make("TOPK", 2, 500.0, 5.0, 5.4, "G1", quantity: 9e6),
            Make("LOWK", 2, 700.0, 8.0, 8.4, "G1", quantity: 1e5),
        };
        var p = Params("LOWK") with { MaxPeptidesPerProtein = 1 };

        var result = Scheduler.Run(cands, p);

        Assert.Equal(new HashSet<string> { "LOWK2" }, Scheduled(result, cands));
        Assert.Equal(ForcedPeptideStatus.Scheduled, result.ForcedOutcomes[0].Status);
    }

    [Fact]
    public void Scheduler_ForcedPeptide_ClaimsBudgetBeforeTheFill()
    {
        // Budget 1 and two co-eluting proteins: the higher-priority protein
        // would normally take the only slot.
        var cands = new[]
        {
            Make("BIGK", 2, 500.0, 5.0, 5.4, "G1", quantity: 9e6),
            Make("SMALLK", 2, 700.0, 5.0, 5.4, "G2", quantity: 1e5),
        };

        var result = Scheduler.Run(cands, Params("SMALLK", budget: 1));

        Assert.Equal(new HashSet<string> { "SMALLK2" }, Scheduled(result, cands));
    }

    [Fact]
    public void Scheduler_ForcedPeptide_FallsBackToAnotherChargeWhenBestFailsTheRules()
    {
        // The 2+ (best-ranked) sits under a dominating co-isolated background
        // and fails the share floor; the 3+ is in clean space.
        var cands = new[]
        {
            Make("PEPK", 2, 500.0, 5.0, 5.4, "G1", quantity: 1e6),
            Make("PEPK", 3, 333.7, 5.0, 5.4, "G1", quantity: 5e5, fragments: new[] { 110.0, 210.0, 310.0, 410.0 }),
            Make("BGK", 2, 500.2, 5.0, 5.4, "BG", quantity: 1e9, fragments: new[] { 610.0, 620.0, 630.0, 640.0 }),
        };
        var p = Params("PEPK", share: 0.05) with
        {
            TargetProteins = new HashSet<string> { "G1" },
            TargetMode = TargetListMode.Exclusive,
        };

        var result = Scheduler.Run(cands, p);

        var o = result.ForcedOutcomes[0];
        Assert.Equal(ForcedPeptideStatus.Scheduled, o.Status);
        Assert.Equal(3, cands[o.CandidateIndex].PrecursorCharge);
        Assert.True(o.UsedFallbackCharge);
    }

    [Fact]
    public void Scheduler_ForcedPeptide_ThatFailsEveryRuleIsReportedRejected()
    {
        var cands = new[]
        {
            Make("PEPK", 2, 500.0, 5.0, 5.4, "G1", quantity: 1e6),
            Make("BGK", 2, 500.2, 5.0, 5.4, "BG", quantity: 1e6, fragments: new[] { 100.1, 620.0, 630.0, 640.0 }),
        };
        var p = Params("PEPK") with
        {
            TargetProteins = new HashSet<string> { "G1" },
            TargetMode = TargetListMode.Exclusive,
        };

        var result = Scheduler.Run(cands, p);

        Assert.Empty(result.ScheduledIndices);
        Assert.Equal(ForcedPeptideStatus.Rejected, result.ForcedOutcomes[0].Status);
    }

    [Fact]
    public void Scheduler_ForcedPeptide_NotInThePoolIsReported()
    {
        var cands = new[] { Make("PEPK", 2, 500.0, 5.0, 5.4, "G1") };

        var result = Scheduler.Run(cands, Params("MISSINGK"));

        Assert.Equal(ForcedPeptideStatus.NotInPool, result.ForcedOutcomes[0].Status);
        Assert.Single(result.ScheduledIndices);
    }

    [Fact]
    public void Scheduler_ForcedWindow_AcceptsFillThatPassesTheRules()
    {
        // Forced windows are ordinary windows: "other" multiplexes into the
        // forced peptide's window exactly as it would without forcing.
        var cands = new[]
        {
            Make("FORCEDK", 2, 500.0, 5.0, 5.4, "G1"),
            Make("OTHERK", 2, 501.0, 5.1, 5.5, "G2", fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
        };

        var result = Scheduler.Run(cands, Params("FORCEDK"));

        Assert.Equal(2, result.ScheduledIndices.Length);
        Assert.Single(result.Slots);
    }

    [Fact]
    public void Scheduler_ForcedPeptides_ArePlacedTheNormalWay()
    {
        // Two compatible forced peptides share a window, as the normal
        // placement would do, rather than each opening its own.
        var cands = new[]
        {
            Make("AAK", 2, 500.0, 5.0, 5.4, "G1"),
            Make("CCK", 2, 501.0, 5.1, 5.5, "G2", fragments: new[] { 150.0, 250.0, 350.0, 450.0 }),
        };

        var result = Scheduler.Run(cands, Params("AAK\nCCK"));

        Assert.Single(result.Slots);
        Assert.All(result.ForcedOutcomes, o => Assert.Equal(ForcedPeptideStatus.Scheduled, o.Status));
    }

    [Fact]
    public void Scheduler_ForcedWindow_RejectsFillThatWouldBreakTheForcedPeptidesRules()
    {
        // "other" would put its 100.1 fragment on the forced peptide's 100.0
        // fragment, so it cannot join the forced window and gets its own.
        var cands = new[]
        {
            Make("FORCEDK", 2, 500.0, 5.0, 5.4, "G1"),
            Make("OTHERK", 2, 501.0, 5.1, 5.5, "G2", fragments: new[] { 100.1, 250.0, 350.0, 450.0 }),
        };

        var result = Scheduler.Run(cands, Params("FORCEDK"));

        Assert.Contains(result.ForcedOutcomes, o => o.Status == ForcedPeptideStatus.Scheduled);
        Assert.DoesNotContain(result.Slots, sl => sl.MemberIndices.Count > 1);
    }

    [Fact]
    public void Scheduler_ForcedPeptide_IgnoresAnExclusiveTargetList()
    {
        var cands = new[]
        {
            Make("INLISTK", 2, 500.0, 5.0, 5.4, "G1"),
            Make("OFFLISTK", 2, 700.0, 8.0, 8.4, "G2"),
        };
        var p = Params("OFFLISTK") with
        {
            TargetProteins = new HashSet<string> { "G1" },
            TargetMode = TargetListMode.Exclusive,
        };

        var result = Scheduler.Run(cands, p);

        Assert.Equal(new HashSet<string> { "INLISTK2", "OFFLISTK2" }, Scheduled(result, cands));
    }

    [Fact]
    public void Scheduler_ForcedPeptide_SurvivesMinPeptidesPerProtein()
    {
        var cands = new[] { Make("LONEK", 2, 500.0, 5.0, 5.4, "G1") };

        var result = Scheduler.Run(cands, Params("LONEK") with { MinPeptidesPerProtein = 2 });

        Assert.Single(result.ScheduledIndices);
    }
}
