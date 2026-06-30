#nullable enable

using System.Globalization;

namespace SkylineCadenza.Core.Ingest;

/// <summary>
/// Reads a predicted peptide-ranking overlay and exposes a per-precursor
/// priority score for Carafe-only scheduling. Two interchangeable source
/// formats are auto-detected from the header:
///
/// <list type="bullet">
/// <item>
/// <b>Bromo</b> <c>peptide_rankings.tsv</c>:
/// <c>protein, peptide (SEQ|charge), mean_pred_score, wins, rank</c>. The
/// score is <c>mean_pred_score</c> (higher = better). Peptides are unmodified,
/// so only the stripped-sequence key is populated.
/// </item>
/// <item>
/// <b>Carafe</b> <c>*_detectability_pred.tsv</c>:
/// <c>sequence, charge, ModifiedPeptide, detectability_pred, ...</c>. The score
/// is <c>detectability_pred</c> (higher = better). Both the modified-sequence
/// key (Carafe <c>_..._</c> form) and the stripped-sequence key are populated.
/// </item>
/// </list>
///
/// The overlay carries no fragments, m/z, or RT, so it is matched onto an
/// existing Carafe-library candidate pool by (peptide, charge) in
/// <see cref="CandidateBuilder.BuildFromCarafe"/>.
/// </summary>
/// <remarks>
/// Same hand-rolled, header-indexed TSV style as
/// <see cref="CarafeLibraryReader"/> - no CSV library, span-based splitting.
/// Scores are NOT normalized: within-protein ordering only needs relative
/// order, and both sources are higher-is-better.
/// </remarks>
public static class PredictedRankingReader
{
    private const string DetectabilitySuffix = "_detectability_pred.tsv";

    // Carafe detectability columns.
    private const string Col_Sequence = "sequence";
    private const string Col_Charge = "charge";
    private const string Col_ModifiedPeptide = "ModifiedPeptide";
    private const string Col_Detectability = "detectability_pred";

    // Bromo peptide-ranking columns.
    private const string Col_Peptide = "peptide";
    private const string Col_MeanPredScore = "mean_pred_score";
    private const string Col_Rank = "rank";

    /// <summary>
    /// A resolved predicted-ranking overlay. Lookups try the modified-sequence
    /// key first, then fall back to the stripped-sequence key.
    /// </summary>
    public sealed class Result
    {
        /// <summary>Keyed by Carafe-form modified peptide + charge.</summary>
        public Dictionary<FragmentKey, double> ByModifiedKey { get; } = new();

        /// <summary>Keyed by stripped sequence + charge.</summary>
        public Dictionary<(string Sequence, int Charge), double> ByStrippedKey { get; } = new();

        /// <summary>Distinct source files merged into this overlay.</summary>
        public List<string> SourceFiles { get; } = new();

        /// <summary>Total distinct (peptide, charge) keys resolved.</summary>
        public int KeyCount => Math.Max(ByModifiedKey.Count, ByStrippedKey.Count);

        /// <summary>
        /// Resolve the predicted score for a precursor, preferring the
        /// modified-sequence match. Returns null when neither key is present.
        /// </summary>
        public double? Lookup(string modifiedPeptide, string strippedPeptide, int charge)
        {
            if (ByModifiedKey.TryGetValue(new FragmentKey(modifiedPeptide, charge), out var m))
                return m;
            if (ByStrippedKey.TryGetValue((strippedPeptide, charge), out var s))
                return s;
            return null;
        }
    }

    /// <summary>
    /// Load an overlay from a single file or folder.
    /// <list type="bullet">
    /// <item>A <b>directory</b> merges every <c>*_detectability_pred.tsv</c> in it.</item>
    /// <item>A <c>*_detectability_pred.tsv</c> <b>file</b> also merges its
    /// siblings in the same folder (so pointing at any one detectability file
    /// loads the full set).</item>
    /// <item>Any other <b>file</b> (e.g. a Bromo <c>peptide_rankings.tsv</c>) is
    /// read on its own.</item>
    /// </list>
    /// On key collision across merged files, the higher score wins.
    /// </summary>
    public static Result LoadFromPath(string path)
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Predicted-ranking path is empty.", nameof(path));
        return Load(ResolveFiles(path, expandSiblings: true), path);
    }

    /// <summary>
    /// Load an overlay from several explicitly chosen files and/or folders,
    /// merging them all into one result (higher score wins on key collision).
    /// Each entry may be a file (Bromo or a single detectability shard) or a
    /// directory (every <c>*_detectability_pred.tsv</c> in it). Use this when
    /// the user multi-selects detectability shards.
    /// </summary>
    /// <remarks>
    /// Unlike <see cref="LoadFromPath"/>, explicitly listed detectability
    /// files are taken literally and do NOT auto-expand to their folder
    /// siblings, so the merge contains exactly what the user selected. (A
    /// single-element call delegates to <see cref="LoadFromPath"/> to preserve
    /// the one-file convenience.) Duplicate files across entries are read once.
    /// </remarks>
    public static Result LoadFromPaths(IEnumerable<string> paths)
    {
        var list = (paths ?? Enumerable.Empty<string>())
            .Where(p => !string.IsNullOrWhiteSpace(p))
            .Select(p => p.Trim())
            .ToList();
        if (list.Count == 0)
            throw new ArgumentException("No predicted-ranking paths supplied.", nameof(paths));
        if (list.Count == 1)
            return LoadFromPath(list[0]);

        var files = new List<string>();
        foreach (var p in list)
            files.AddRange(ResolveFiles(p, expandSiblings: false));
        return Load(files, string.Join(", ", list));
    }

    private static Result Load(List<string> files, string origin)
    {
        var result = new Result();
        var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (var file in files)
        {
            string full = Path.GetFullPath(file);
            if (!seen.Add(full)) continue; // read each file once
            ReadInto(file, result);
            result.SourceFiles.Add(file);
        }
        if (result.SourceFiles.Count == 0)
            throw new FileNotFoundException(
                $"No predicted-ranking files found at '{origin}'.", origin);
        return result;
    }

    private static List<string> ResolveFiles(string path, bool expandSiblings)
    {
        if (Directory.Exists(path))
            return Directory.GetFiles(path, "*" + DetectabilitySuffix)
                .OrderBy(f => f, StringComparer.OrdinalIgnoreCase).ToList();

        if (!File.Exists(path))
            return new List<string>();

        // A lone detectability file pulls in its siblings; an explicitly listed
        // one (expandSiblings == false) and any other file stand alone.
        if (expandSiblings
            && path.EndsWith(DetectabilitySuffix, StringComparison.OrdinalIgnoreCase))
        {
            string? dir = Path.GetDirectoryName(path);
            if (!string.IsNullOrEmpty(dir))
                return Directory.GetFiles(dir, "*" + DetectabilitySuffix)
                    .OrderBy(f => f, StringComparer.OrdinalIgnoreCase).ToList();
        }
        return new List<string> { path };
    }

    private static void ReadInto(string file, Result result)
    {
        using var reader = new StreamReader(File.OpenRead(file));
        string? headerLine = reader.ReadLine()
            ?? throw new InvalidDataException($"Predicted-ranking file is empty: {file}");
        var headers = headerLine.Split('\t');

        int idxDetect = IndexOf(headers, Col_Detectability);
        if (idxDetect >= 0)
        {
            ReadDetectability(reader, headers, idxDetect, result);
            return;
        }

        int idxScore = IndexOf(headers, Col_MeanPredScore);
        int idxPeptide = IndexOf(headers, Col_Peptide);
        if (idxScore >= 0 && idxPeptide >= 0)
        {
            ReadBromo(reader, idxPeptide, idxScore, result);
            return;
        }

        throw new InvalidDataException(
            $"Unrecognized predicted-ranking format in '{file}'. Expected a Carafe " +
            $"'{Col_Detectability}' column or Bromo '{Col_MeanPredScore}'/'{Col_Peptide}' " +
            $"columns. Found: {string.Join(", ", headers)}");
    }

    private static void ReadDetectability(
        StreamReader reader, string[] headers, int idxDetect, Result result)
    {
        int idxMod = IndexOf(headers, Col_ModifiedPeptide);
        int idxSeq = IndexOf(headers, Col_Sequence);
        int idxCharge = IndexOf(headers, Col_Charge);
        int maxIdx = new[] { idxDetect, idxMod, idxSeq, idxCharge }.Max();

        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            var fields = line.Split('\t');
            if (fields.Length <= maxIdx) continue;
            if (!int.TryParse(fields[idxCharge], NumberStyles.Integer,
                CultureInfo.InvariantCulture, out int charge)) continue;
            if (!double.TryParse(fields[idxDetect], NumberStyles.Float,
                CultureInfo.InvariantCulture, out double score)) continue;

            if (idxMod >= 0)
            {
                string mod = fields[idxMod];
                if (mod.Length > 0)
                    Keep(result.ByModifiedKey, new FragmentKey(mod, charge), score);
            }
            if (idxSeq >= 0)
            {
                string seq = fields[idxSeq];
                if (seq.Length > 0)
                    Keep(result.ByStrippedKey, (seq, charge), score);
            }
        }
    }

    private static void ReadBromo(
        StreamReader reader, int idxPeptide, int idxScore, Result result)
    {
        int maxIdx = Math.Max(idxPeptide, idxScore);
        string? line;
        while ((line = reader.ReadLine()) != null)
        {
            var fields = line.Split('\t');
            if (fields.Length <= maxIdx) continue;
            if (!double.TryParse(fields[idxScore], NumberStyles.Float,
                CultureInfo.InvariantCulture, out double score)) continue;

            // peptide is "SEQ|charge"; split on the last '|'.
            string peptide = fields[idxPeptide];
            int bar = peptide.LastIndexOf('|');
            if (bar <= 0 || bar == peptide.Length - 1) continue;
            string seq = peptide.Substring(0, bar);
            if (!int.TryParse(peptide.AsSpan(bar + 1), NumberStyles.Integer,
                CultureInfo.InvariantCulture, out int charge)) continue;

            Keep(result.ByStrippedKey, (seq, charge), score);
        }
    }

    private static void Keep<TKey>(Dictionary<TKey, double> map, TKey key, double score)
        where TKey : notnull
    {
        if (!map.TryGetValue(key, out var existing) || score > existing)
            map[key] = score;
    }

    private static int IndexOf(string[] headers, string name)
    {
        for (int i = 0; i < headers.Length; i++)
            if (headers[i] == name) return i;
        return -1;
    }
}
