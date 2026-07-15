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

    private static DiannRow DiannRow(string strippedSeq, int charge, string proteinId) =>
        new()
        {
            PrecursorId = $"{strippedSeq}+{charge}",
            ModifiedSequence = strippedSeq,
            StrippedSequence = strippedSeq,
            PrecursorCharge = charge,
            PrecursorMz = 500.0,
            RtApex = 20.0,
            RtStart = 19.8,
            RtStop = 20.2,
            PrecursorQuantity = 1e6,
            Proteotypic = 1,
            ProteinIds = proteinId,
            Genes = "",
            Run = "run1",
            QValue = 0.001,
            ProteinQValue = 0.001,
            Fragments = System.Array.Empty<FragmentIon>(),
        };

    private static CarafePrecursor Carafe(string strippedSeq, int charge, string accession) =>
        new()
        {
            ModifiedPeptide = $"_{strippedSeq}_",
            StrippedPeptide = strippedSeq,
            PrecursorMz = 500.0,
            PrecursorCharge = charge,
            TrRecalibrated = 20.0,
            ProteinIdRaw = $"sp|{accession}|NAME",
            ProteinAccessions = new[] { accession },
            PredictedIntensity = 1.0,
            Fragments = System.Array.Empty<FragmentIon>(),
        };

    [Fact]
    public void BuildCombined_DedupsPredictedDuplicateOnPeptideAndCharge_KeepsOtherChargeAndNewPeptide()
    {
        // Confirmed by DIA-NN: PEPTIDEK+2.
        var rows = new[] { DiannRow("PEPTIDEK", 2, "P1") };
        var carafe = new[]
        {
            Carafe("PEPTIDEK", 2, "P1"), // exact duplicate of the confirmed target -> dropped
            Carafe("PEPTIDEK", 3, "P1"), // same peptide, different charge -> kept as predicted
            Carafe("OTHERPEP", 2, "P2"), // brand-new peptide -> kept as predicted
        };
        var parsimony = new Dictionary<string, (string Group, string Type)>
        {
            ["PEPTIDEK"] = ("P1", "unique"),
            ["OTHERPEP"] = ("P2", "unique"),
        };

        var result = CandidateBuilder.BuildCombined(
            rows, new Dictionary<FragmentKey, FragmentIon[]>(), carafe, parsimony);

        // 1 confirmed + 2 predicted survivors.
        Assert.Equal(3, result.Count);

        var confirmed = result.Where(c => c.Source == CandidateSource.Experimental).ToList();
        Assert.Equal("PEPTIDEK", Assert.Single(confirmed).StrippedSequence);
        Assert.Equal(2, confirmed[0].PrecursorCharge);

        var predicted = result.Where(c => c.Source == CandidateSource.Predicted)
            .Select(c => (c.StrippedSequence, c.PrecursorCharge)).ToHashSet();
        Assert.DoesNotContain(("PEPTIDEK", 2), predicted); // deduped
        Assert.Contains(("PEPTIDEK", 3), predicted);       // different charge kept
        Assert.Contains(("OTHERPEP", 2), predicted);       // new peptide kept
    }
}
