using SkylineCadenza.Core.Ingest;
using Xunit;

namespace SkylineCadenza.Tests;

public class PredictedRankingReaderTests
{
    private static string WriteTemp(string fileName, string content)
    {
        string dir = Path.Combine(Path.GetTempPath(), "cadenza_pr_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        string path = Path.Combine(dir, fileName);
        File.WriteAllText(path, content);
        return path;
    }

    [Fact]
    public void Bromo_Format_KeysByStrippedSequenceAndCharge()
    {
        string path = WriteTemp("peptide_rankings.tsv",
            "protein\tpeptide\tmean_pred_score\twins\trank\n" +
            "P1\tALIAVSEPEAR|2\t0.95\t538\t1\n" +
            "P1\tVTSQVVDR|3\t0.88\t522\t2\n");

        var r = PredictedRankingReader.LoadFromPath(path);

        Assert.Equal(0.95, r.Lookup("anything", "ALIAVSEPEAR", 2));
        Assert.Equal(0.88, r.Lookup("anything", "VTSQVVDR", 3));
        // Charge must match.
        Assert.Null(r.Lookup("anything", "ALIAVSEPEAR", 3));
        // Bromo populates the stripped key only.
        Assert.Empty(r.ByModifiedKey);
    }

    [Fact]
    public void Detectability_Format_KeysByModifiedAndStripped()
    {
        string path = WriteTemp("1_detectability_pred.tsv",
            "sequence\tmods\tmod_sites\tcharge\tModifiedPeptide\tdetectability_pred\n" +
            "NWEYNPLYPREWPR\t\t\t3\t_NWEYNPLYPREWPR_\t0.5993\n" +
            "FQCDNHRCIPEQWK\tCarbamidomethyl@C\t3\t4\t_FQC[UniMod:4]DNHRCIPEQWK_\t0.2259\n");

        var r = PredictedRankingReader.LoadFromPath(path);

        // Modified key preferred.
        Assert.Equal(0.5993, r.Lookup("_NWEYNPLYPREWPR_", "NWEYNPLYPREWPR", 3));
        Assert.Equal(0.2259, r.Lookup("_FQC[UniMod:4]DNHRCIPEQWK_", "FQCDNHRCIPEQWK", 4));
        // Stripped fallback still resolves when the modified form is unknown.
        Assert.Equal(0.5993, r.Lookup("unknown-mod-form", "NWEYNPLYPREWPR", 3));
    }

    [Fact]
    public void Folder_MergesSiblingDetectabilityFiles_KeepingMaxOnCollision()
    {
        string p0 = WriteTemp("0_detectability_pred.tsv",
            "sequence\tcharge\tModifiedPeptide\tdetectability_pred\n" +
            "PEPTIDEK\t2\t_PEPTIDEK_\t0.30\n");
        // Second file in the SAME folder, with a higher score for the same key.
        string dir = Path.GetDirectoryName(p0)!;
        File.WriteAllText(Path.Combine(dir, "1_detectability_pred.tsv"),
            "sequence\tcharge\tModifiedPeptide\tdetectability_pred\n" +
            "PEPTIDEK\t2\t_PEPTIDEK_\t0.70\n" +
            "OTHERR\t2\t_OTHERR_\t0.40\n");

        // Pointing at one detectability file pulls in its siblings.
        var r = PredictedRankingReader.LoadFromPath(p0);

        Assert.Equal(2, r.SourceFiles.Count);
        Assert.Equal(0.70, r.Lookup("_PEPTIDEK_", "PEPTIDEK", 2)); // max wins
        Assert.Equal(0.40, r.Lookup("_OTHERR_", "OTHERR", 2));
    }

    [Fact]
    public void Directory_Path_LoadsAllDetectabilityFiles()
    {
        string p0 = WriteTemp("0_detectability_pred.tsv",
            "sequence\tcharge\tModifiedPeptide\tdetectability_pred\n" +
            "AAAAK\t2\t_AAAAK_\t0.11\n");
        string dir = Path.GetDirectoryName(p0)!;

        var r = PredictedRankingReader.LoadFromPath(dir);

        Assert.Equal(0.11, r.Lookup("_AAAAK_", "AAAAK", 2));
    }

    [Fact]
    public void MultipleFiles_MergeExactlySelected_NoSiblingExpansion()
    {
        // Three shards in one folder; the user explicitly selects only two.
        string s0 = WriteTemp("0_detectability_pred.tsv",
            "sequence\tcharge\tModifiedPeptide\tdetectability_pred\n" +
            "AAA\t2\t_AAA_\t0.10\n");
        string dir = Path.GetDirectoryName(s0)!;
        string s1 = Path.Combine(dir, "1_detectability_pred.tsv");
        File.WriteAllText(s1,
            "sequence\tcharge\tModifiedPeptide\tdetectability_pred\n" +
            "BBB\t2\t_BBB_\t0.20\n");
        string s2 = Path.Combine(dir, "2_detectability_pred.tsv");
        File.WriteAllText(s2,
            "sequence\tcharge\tModifiedPeptide\tdetectability_pred\n" +
            "CCC\t2\t_CCC_\t0.30\n");

        var r = PredictedRankingReader.LoadFromPaths(new[] { s0, s1 });

        Assert.Equal(2, r.SourceFiles.Count);
        Assert.Equal(0.10, r.Lookup("_AAA_", "AAA", 2));
        Assert.Equal(0.20, r.Lookup("_BBB_", "BBB", 2));
        // The unselected sibling must NOT be pulled in.
        Assert.Null(r.Lookup("_CCC_", "CCC", 2));
    }

    [Fact]
    public void MultipleFiles_DedupRepeatedPathAndKeepMax()
    {
        string a = WriteTemp("bromo_a.tsv",
            "protein\tpeptide\tmean_pred_score\twins\trank\n" +
            "P1\tPEPK|2\t0.40\t1\t1\n");
        string b = WriteTemp("bromo_b.tsv",
            "protein\tpeptide\tmean_pred_score\twins\trank\n" +
            "P1\tPEPK|2\t0.80\t1\t1\n");

        // 'a' listed twice: it must be read only once, and the higher score
        // across the distinct files wins.
        var r = PredictedRankingReader.LoadFromPaths(new[] { a, a, b });

        Assert.Equal(2, r.SourceFiles.Count);
        Assert.Equal(0.80, r.Lookup("x", "PEPK", 2));
    }

    [Fact]
    public void UnknownHeader_Throws()
    {
        string path = WriteTemp("mystery.tsv", "foo\tbar\nx\ty\n");
        Assert.Throws<InvalidDataException>(() => PredictedRankingReader.LoadFromPath(path));
    }
}
