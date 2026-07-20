using SkylineCadenza.Core.Output;
using SkylineCadenza.Core.Scheduling;
using Xunit;

namespace SkylineCadenza.Tests;

public class PeptideTransitionListBuilderTests
{
    private static Candidate Make(
        string seq, string modSeq, int z, double mz,
        double rtApex, string protein,
        params FragmentIon[] frags) =>
        new Candidate
        {
            PrecursorId = $"{seq}+{z}",
            StrippedSequence = seq,
            ModifiedSequence = modSeq,
            PrecursorCharge = z,
            PrecursorMz = mz,
            RtStart = rtApex - 0.25,
            RtStop = rtApex + 0.25,
            RtApex = rtApex,
            PrecursorQuantity = 1e6,
            QValue = 0.001,
            ProteinQValue = 0.001,
            Proteotypic = 1,
            ProteinGroup = protein,
            PeptideType = "unique",
            Fragments = frags,
            Top4Fragments = Candidate.DeriveTopMz(frags, 4),
            Run = "test",
        };

    private static ScheduleResult OneSlot(IReadOnlyList<Candidate> cands)
    {
        var slot = new Slot
        {
            Id = 0,
            MzMin = cands.Min(c => c.PrecursorMz),
            MzMax = cands.Max(c => c.PrecursorMz),
            RtStart = cands.Min(c => c.RtStart),
            RtStop = cands.Max(c => c.RtStop),
            CoStart = cands.Min(c => c.RtStart),
            CoStop = cands.Max(c => c.RtStop),
            Fragments = Array.Empty<double>(),
        };
        for (int i = 0; i < cands.Count; i++) slot.MemberIndices.Add(i);

        return new ScheduleResult
        {
            ScheduledIndices = Enumerable.Range(0, cands.Count).ToArray(),
            ScheduledSlotIds = new int[cands.Count], // all zeros = slot 0
            Slots = new[] { slot },
            RtGrid = Array.Empty<double>(),
            SlotCountCurve = Array.Empty<int>(),
            ProteinGroupsCovered = cands.Select(c => c.ProteinGroup).Distinct().Count(),
        };
    }

    [Fact]
    public void Build_EmitsRealProductCharge_PerFragment()
    {
        // Mix of +1 and +2 fragments - the previous hardcoded +1 builder
        // would have emitted +1 for both and Skyline would reject the +2
        // m/z values as "no matching product ion".
        var cands = new[]
        {
            Make("ASFNHFDK", "ASFNHFDK", 2, 471.71, 5.0, "P1",
                new FragmentIon(404.1928, 100.0, 2),  // y6+2
                new FragmentIon(447.7089, 80.0, 2)),  // y7+2
        };
        var schedule = OneSlot(cands);
        var csv = PeptideTransitionListBuilder.Build(cands, schedule);
        var lines = csv.Trim().Split('\n');

        // Header + 3 precursor isotope rows + 2 fragment rows.
        Assert.Equal(1 + PeptideTransitionListBuilder.PrecursorIsotopesPerPeptide + 2,
            lines.Length);
        Assert.StartsWith("Protein Name,Peptide Modified Sequence,", lines[0]);

        // Fragment rows follow the precursor isotope rows; check their
        // Product Charge to confirm the +2 fragments survive.
        var fragRow1 = lines[1 + PeptideTransitionListBuilder.PrecursorIsotopesPerPeptide]
            .Split(',');
        var fragRow2 = lines[2 + PeptideTransitionListBuilder.PrecursorIsotopesPerPeptide]
            .Split(',');
        Assert.Equal("404.1928", fragRow1[4]); // Product m/z
        Assert.Equal("2", fragRow1[5]);        // Product Charge
        Assert.Equal("447.7089", fragRow2[4]);
        Assert.Equal("2", fragRow2[5]);
    }

    [Fact]
    public void Build_EmitsPrecursorIsotopeRows_M0M1M2_PerPeptide()
    {
        // Skyline detects a precursor transition when Product Charge ==
        // Precursor Charge and Product m/z is at the precursor's
        // monoisotope spacing. M+0 = precursor m/z; M+1 = +1.003355 / z;
        // M+2 = +2 * 1.003355 / z.
        const int charge = 2;
        const double precursorMz = 500.0000;
        var cands = new[]
        {
            Make("PEPTIDE", "PEPTIDE", charge, precursorMz, 5.0, "P1",
                new FragmentIon(200.0, 100.0, 1)),
        };
        var schedule = OneSlot(cands);
        var lines = PeptideTransitionListBuilder.Build(cands, schedule)
            .Trim().Split('\n');

        // Header + 3 precursor isotope rows + 1 fragment row.
        Assert.Equal(1 + PeptideTransitionListBuilder.PrecursorIsotopesPerPeptide + 1,
            lines.Length);

        const double neutronMass = 1.003355;
        for (int isotope = 0; isotope < PeptideTransitionListBuilder.PrecursorIsotopesPerPeptide;
             isotope++)
        {
            var fields = lines[1 + isotope].Split(',');
            double expectedMz = precursorMz + neutronMass * isotope / charge;
            Assert.Equal(expectedMz, double.Parse(fields[4],
                System.Globalization.CultureInfo.InvariantCulture), precision: 4);
            Assert.Equal(charge.ToString(), fields[5]); // Product Charge == Precursor Charge
        }
    }

    [Fact]
    public void Build_LimitsToSixFragmentsPerPeptide()
    {
        var frags = Enumerable.Range(1, 10)
            .Select(i => new FragmentIon(100.0 + i, 20.0 - i, i % 2 == 0 ? 2 : 1))
            .ToArray();
        var cands = new[] { Make("MANYFRAG", "MANYFRAG", 2, 400.0, 5.0, "P1", frags) };
        var schedule = OneSlot(cands);
        var csv = PeptideTransitionListBuilder.Build(cands, schedule);

        int rowCount = csv.Trim().Split('\n').Length - 1; // exclude header
        // Precursor isotope rows are emitted in addition to the
        // fragment cap.
        Assert.Equal(
            PeptideTransitionListBuilder.PrecursorIsotopesPerPeptide
            + PeptideTransitionListBuilder.FragmentsPerPeptide,
            rowCount);
    }

    [Fact]
    public void Build_SkipsHeavyLabeledPrecursors_KeepsLight()
    {
        // Heavy isotope-labeled precursors (e.g. PRTC standards) carry the
        // document's explicit isotope modification, which Skyline's
        // headless transition-list importer rejects and which - because
        // the row is used to validate the precursor-m/z column - aborts
        // the entire import. The builder omits them; they reach the assay
        // via the registered library instead.
        var light = Make("PEPTIDER", "PEPTIDER", 2, 400.0, 5.0, "P1",
            new FragmentIon(200.0, 100.0, 1));
        var heavy = new Candidate
        {
            PrecursorId = "PRTCPEP+2",
            StrippedSequence = "PRTCPEP",
            ModifiedSequence = "PRTCPEP",
            PrecursorCharge = 2,
            PrecursorMz = 500.0,
            RtStart = 4.75,
            RtStop = 5.25,
            RtApex = 5.0,
            PrecursorQuantity = 1e6,
            QValue = 0.001,
            ProteinQValue = 0.001,
            Proteotypic = 1,
            ProteinGroup = "PRTC",
            PeptideType = "unique",
            Fragments = new[] { new FragmentIon(250.0, 100.0, 1) },
            Top4Fragments = Candidate.DeriveTopMz(new[] { new FragmentIon(250.0, 100.0, 1) }, 4),
            Run = "test",
            IsotopeLabelType = "heavy",
        };
        var cands = new[] { light, heavy };
        var schedule = OneSlot(cands);

        var csv = PeptideTransitionListBuilder.Build(cands, schedule);

        // Only the light peptide's rows survive; the heavy one contributes none.
        Assert.Contains("PEPTIDER", csv);
        Assert.DoesNotContain("PRTCPEP", csv);
    }

    [Fact]
    public void Build_FragmentlessPeptide_StillEmitsPrecursorIsotopeRows()
    {
        var cands = new[] { Make("NOFRAGS", "NOFRAGS", 2, 400.0, 5.0, "P1") };
        var schedule = OneSlot(cands);
        var csv = PeptideTransitionListBuilder.Build(cands, schedule);
        var lines = csv.Trim().Split('\n');

        // Header + 3 precursor isotope rows + 1 precursor-only row
        // (the fragmentless fallback).
        Assert.Equal(1 + PeptideTransitionListBuilder.PrecursorIsotopesPerPeptide + 1,
            lines.Length);
        var fallback = lines[^1].Split(',');
        Assert.Equal(string.Empty, fallback[4]); // Product m/z blank
        Assert.Equal(string.Empty, fallback[5]); // Product Charge blank
    }
}
