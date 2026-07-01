using SkylineCadenza.Core.Ingest;
using SkylineCadenza.Core.Scheduling;
using Xunit;

namespace SkylineCadenza.Tests;

public class CandidateBuilderTests
{
    private static List<CarafePrecursor> OnePrecursor() => new()
    {
        new CarafePrecursor
        {
            ModifiedPeptide = "_PEPTIDEK_",
            StrippedPeptide = "PEPTIDEK",
            PrecursorMz = 500.0,
            PrecursorCharge = 2,
            TrRecalibrated = 20.0,
            ProteinIdRaw = "sp|P1|NAME",
            ProteinAccessions = new[] { "P1" },
            PredictedIntensity = 1.0,
            Fragments = System.Array.Empty<FragmentIon>(),
        },
    };

    private static Dictionary<string, (string Group, string Type)> Parsimony() =>
        new() { ["PEPTIDEK"] = ("P1", "unique") };

    private static PredictedRankingReader.Result Overlay(double score)
    {
        var r = new PredictedRankingReader.Result();
        r.ByStrippedKey[("PEPTIDEK", 2)] = score;
        return r;
    }

    [Fact]
    public void BuildFromCarafe_NoOverlays_PredictedPriorityIsNull()
    {
        var c = CandidateBuilder.BuildFromCarafe(OnePrecursor(), Parsimony());
        Assert.Null(Assert.Single(c).PredictedPriority);
    }

    [Fact]
    public void BuildFromCarafe_DetectabilityOnly_UsedAsRankingFallback()
    {
        // No Bromo overlay; detectability supplied with no cutoff -> it should
        // populate PredictedPriority (ranking fallback), not just filter.
        var c = CandidateBuilder.BuildFromCarafe(
            OnePrecursor(), Parsimony(), predicted: null,
            detectability: Overlay(0.87), minDetectability: 0.0);
        Assert.Equal(0.87, Assert.Single(c).PredictedPriority);
    }

    [Fact]
    public void BuildFromCarafe_BromoTakesPrecedenceOverDetectability()
    {
        var c = CandidateBuilder.BuildFromCarafe(
            OnePrecursor(), Parsimony(), predicted: Overlay(0.50),
            detectability: Overlay(0.87), minDetectability: 0.0);
        Assert.Equal(0.50, Assert.Single(c).PredictedPriority); // Bromo wins
    }

    [Fact]
    public void BuildFromCarafe_DetectabilityCutoff_DropsBelowThreshold()
    {
        var c = CandidateBuilder.BuildFromCarafe(
            OnePrecursor(), Parsimony(), predicted: null,
            detectability: Overlay(0.30), minDetectability: 0.5);
        Assert.Empty(c);
    }

    [Fact]
    public void BuildFromCarafe_DetectabilityCutoff_KeepsUnknown()
    {
        // Peptide not in the detectability overlay -> unknown -> kept.
        var c = CandidateBuilder.BuildFromCarafe(
            OnePrecursor(), Parsimony(), predicted: null,
            detectability: new PredictedRankingReader.Result(), minDetectability: 0.5);
        Assert.Single(c);
        Assert.Null(c[0].PredictedPriority);
    }
}
