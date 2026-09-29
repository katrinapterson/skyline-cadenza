using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;
using SkylineCadenza.Core.Output;

namespace SkylineCadenza.Core.Scheduling;

/// <summary>
/// One entry of the forced-peptide list: a peptide the user wants in the
/// schedule regardless of its ranking, e.g. to carry the same validation
/// peptides across several schedules.
/// </summary>
/// <param name="Text">The entry as the user wrote it (for reporting).</param>
/// <param name="Key">
/// Matching key: the stripped sequence, or the normalized modified sequence
/// when <paramref name="Modified"/> is true (see <see cref="ForcedPeptides.ModifiedKey"/>).
/// </param>
/// <param name="Modified">True when the entry carries modifications and must match them.</param>
/// <param name="Charge">A specific precursor charge, or null for "any charge, best first".</param>
public readonly record struct ForcedPeptideSpec(string Text, string Key, bool Modified, int? Charge);

/// <summary>What happened to one <see cref="ForcedPeptideSpec"/>.</summary>
public enum ForcedPeptideStatus
{
    /// <summary>Scheduled (see <see cref="ForcedPeptideOutcome.CandidateIndex"/>).</summary>
    Scheduled,

    /// <summary>No candidate in the pool matches: no identification, so no RT to schedule at.</summary>
    NotInPool,

    /// <summary>Matching precursors exist, but every one failed the scheduling rules.</summary>
    Rejected,
}

/// <param name="Spec">The list entry.</param>
/// <param name="Status">Outcome.</param>
/// <param name="CandidateIndex">The scheduled precursor, or -1.</param>
/// <param name="UsedFallbackCharge">
/// True when the best-ranked charge state could not be placed and another
/// charge state of the same peptide was scheduled instead.
/// </param>
public readonly record struct ForcedPeptideOutcome(
    ForcedPeptideSpec Spec, ForcedPeptideStatus Status, int CandidateIndex, bool UsedFallbackCharge);

/// <summary>Parsing and matching for the forced-peptide list.</summary>
public static partial class ForcedPeptides
{
    /// <summary>
    /// Parse one entry per line. Accepted forms: <c>PEPTIDEK</c> (any charge),
    /// <c>PEPTIDEK+2</c>, <c>PEPTIDEK/2</c>, <c>PEPTIDEK++</c>, or a charge in a
    /// second comma / tab / space separated column. Modified sequences in
    /// Skyline (<c>C[+57.021464]</c>, <c>C[UniMod:4]</c>), DIA-NN
    /// (<c>C(UniMod:4)</c>) or Carafe (<c>_C[UniMod:4]_</c>) syntax must match
    /// the modifications exactly; an unmodified entry matches every modified
    /// form of that sequence. Blank lines and lines starting with <c>#</c> are
    /// ignored, and so is anything that is not a peptide (a header row, say),
    /// which is returned in <paramref name="invalid"/>. Duplicates are dropped.
    /// </summary>
    public static List<ForcedPeptideSpec> Parse(string? text, out List<string> invalid)
    {
        invalid = new List<string>();
        var specs = new List<ForcedPeptideSpec>();
        if (string.IsNullOrWhiteSpace(text)) return specs;
        var seen = new HashSet<(string, bool, int?)>();
        foreach (string raw in text.Split('\n'))
        {
            string line = raw.Trim().Trim('"');
            if (line.Length == 0 || line.StartsWith('#')) continue;

            string seq = line;
            int? charge = null;
            var cols = line.Split(new[] { ',', '\t', ' ', ';' }, StringSplitOptions.RemoveEmptyEntries);
            if (cols.Length >= 2 && int.TryParse(cols[1], NumberStyles.None, CultureInfo.InvariantCulture, out int zc))
            {
                seq = cols[0];
                charge = zc;
            }
            else
            {
                seq = cols[0];
                var m = ChargeSuffix().Match(seq);
                if (m.Success)
                {
                    charge = m.Groups["z"].Success
                        ? int.Parse(m.Groups["z"].Value, CultureInfo.InvariantCulture)
                        : m.Groups["plus"].Value.Length;
                    seq = seq[..m.Index];
                }
            }

            string stripped = StrippedKey(seq);
            if (stripped.Length == 0 || charge is <= 0 || !IsPeptide(seq))
            {
                invalid.Add(raw.Trim());
                continue;
            }
            bool modified = seq.IndexOfAny(new[] { '[', '(' }) >= 0;
            string key = modified ? ModifiedKey(seq) : stripped;
            if (seen.Add((key, modified, charge)))
                specs.Add(new ForcedPeptideSpec(line, key, modified, charge));
        }
        return specs;
    }

    /// <summary>True when <paramref name="c"/> is a precursor of <paramref name="spec"/>.</summary>
    public static bool Matches(ForcedPeptideSpec spec, Candidate c)
    {
        if (spec.Charge is int z && c.PrecursorCharge != z) return false;
        return spec.Modified
            ? ModifiedKey(c.ModifiedSequence) == spec.Key
            : string.Equals(c.StrippedSequence, spec.Key, StringComparison.OrdinalIgnoreCase);
    }

    /// <summary>Residues only, upper-cased, with every bracketed modification removed.</summary>
    public static string StrippedKey(string seq)
    {
        var sb = new StringBuilder(seq.Length);
        int depth = 0;
        foreach (char ch in seq)
        {
            if (ch is '[' or '(') { depth++; continue; }
            if (ch is ']' or ')') { depth = Math.Max(0, depth - 1); continue; }
            if (depth == 0 && char.IsLetter(ch)) sb.Append(char.ToUpperInvariant(ch));
        }
        return sb.ToString();
    }

    /// <summary>
    /// Modified-sequence key: folded to Skyline's mass-delta form (so DIA-NN,
    /// Carafe and Skyline syntax agree) with every mass rounded to two
    /// decimals, so <c>C[+57.02]</c> and <c>C[+57.021464]</c> are the same key.
    /// </summary>
    public static string ModifiedKey(string modSeq)
    {
        string norm = BlibAssayWriter.NormalizeModifiedSequence(modSeq);
        return MassDelta().Replace(norm, m =>
        {
            if (!double.TryParse(m.Groups["v"].Value, NumberStyles.Float, CultureInfo.InvariantCulture, out double v))
                return m.Value;
            return "[" + (v >= 0 ? "+" : "") + v.ToString("0.00", CultureInfo.InvariantCulture) + "]";
        }).ToUpperInvariant();
    }

    // Only residue letters outside brackets, and brackets that close.
    private static bool IsPeptide(string seq)
    {
        int depth = 0;
        foreach (char ch in seq)
        {
            if (ch is '[' or '(') depth++;
            else if (ch is ']' or ')') { if (--depth < 0) return false; }
            else if (depth == 0 && ch != '_' && !(ch >= 'A' && ch <= 'Z')) return false;
        }
        return depth == 0;
    }

    [GeneratedRegex(@"(?:[+/](?<z>\d+)|(?<plus>\++))$")]
    private static partial Regex ChargeSuffix();

    [GeneratedRegex(@"\[(?<v>[+-]?\d+(?:\.\d+)?)\]")]
    private static partial Regex MassDelta();
}
