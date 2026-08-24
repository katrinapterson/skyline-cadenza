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
    public void Build_PrefersMeasuredDiannFragmentsOverCarafeOverlay()
    {
        // Confirmed peptides must schedule on their MEASURED DIA-NN fragments so
        // loading a Carafe library cannot move the empirical schedule. A Carafe
        // overlay is present for this peptide but must be ignored in favour of
        // the measured fragments.
        var row = new DiannRow
        {
            PrecursorId = "PEPTIDEK+2",
            ModifiedSequence = "PEPTIDEK",
            StrippedSequence = "PEPTIDEK",
            PrecursorCharge = 2,
            PrecursorMz = 500.0,
            RtApex = 20.0,
            RtStart = 19.8,
            RtStop = 20.2,
            PrecursorQuantity = 1e6,
            Proteotypic = 1,
            ProteinIds = "P1",
            Genes = "",
            Run = "run1",
            QValue = 0.001,
            ProteinQValue = 0.001,
            Fragments = new[]
            {
                new FragmentIon(200.0, 9, 1),
                new FragmentIon(300.0, 8, 1),
                new FragmentIon(400.0, 7, 1),
            },
        };
        var carafeFrags = new Dictionary<FragmentKey, FragmentIon[]>
        {
            [new FragmentKey(CarafeKey.FromDiann("PEPTIDEK"), 2)] = new[]
            {
                new FragmentIon(600.0, 9, 1),
                new FragmentIon(700.0, 8, 1),
            },
        };
        var parsimony = new Dictionary<string, (string Group, string Type)>
        {
            ["PEPTIDEK"] = ("P1", "unique"),
        };

        var candidate = Assert.Single(CandidateBuilder.Build(new[] { row }, parsimony, carafeFrags));

        Assert.Contains(200.0, candidate.Top4Fragments);       // measured DIA-NN fragment
        Assert.DoesNotContain(600.0, candidate.Top4Fragments); // Carafe overlay ignored
    }

    [Fact]
    public void Build_FallsBackToCarafeFragmentsWhenDiannHasNone()
    {
        // When DIA-NN carries no measured fragments, the Carafe overlay is still
        // used so the interference check has something to work with.
        var row = DiannRow("PEPTIDEK", 2, "P1"); // helper sets Fragments = empty
        var carafeFrags = new Dictionary<FragmentKey, FragmentIon[]>
        {
            [new FragmentKey(CarafeKey.FromDiann("PEPTIDEK"), 2)] = new[]
            {
                new FragmentIon(600.0, 9, 1),
                new FragmentIon(700.0, 8, 1),
            },
        };
        var parsimony = new Dictionary<string, (string Group, string Type)>
        {
            ["PEPTIDEK"] = ("P1", "unique"),
        };

        var candidate = Assert.Single(CandidateBuilder.Build(new[] { row }, parsimony, carafeFrags));

        Assert.Contains(600.0, candidate.Top4Fragments); // Carafe fallback used
    }

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

    [Fact]
    public void MapPredictedToExperimentalGroups_routes_by_confirmed_accession_and_parsimonises_the_unseen()
    {
        // Confirmed accession map: P1 and P2 both resolve to canonical group P1;
        // P3 is a distinct confirmed group. A Carafe peptide is routed to the
        // group of whichever confirmed protein its accessions hit; peptides on
        // proteins DIA-NN never saw are grouped among themselves.
        var expAccToGroup = new Dictionary<string, string>
        {
            ["P1"] = "P1",
            ["P2"] = "P1",
            ["P3"] = "P3",
        };
        var carafePepProteins = new Dictionary<string, IReadOnlyList<string>>
        {
            ["SEENPEP"] = new[] { "P2" },        // -> confirmed group P1, unique
            ["SPANPEP"] = new[] { "P1", "P3" },  // spans two confirmed groups -> razor
            ["NEWPEP1"] = new[] { "U1" },        // unseen protein
            ["NEWPEP2"] = new[] { "U1" },        // unseen, same protein -> shared canon
        };

        var map = CandidateBuilder.MapPredictedToExperimentalGroups(
            carafePepProteins, expAccToGroup);

        Assert.Equal(("P1", "unique"), map["SEENPEP"]);
        Assert.Equal(("P1", "razor"), map["SPANPEP"]); // min of {P1,P3}
        // Unseen peptides are parsimonised among themselves, not dropped, and do
        // not borrow a confirmed group.
        Assert.Equal("U1", map["NEWPEP1"].Group);
        Assert.Equal("U1", map["NEWPEP2"].Group);
    }
}
