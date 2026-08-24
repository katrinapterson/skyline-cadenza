namespace SkylineCadenza.Core.Parsimony;

/// <summary>
/// PRISM-style protein parsimony. Port of <c>_find_subsumable_fast</c> +
/// <c>build_parsimony_assignment</c> from the targeted-modeling notebook
/// (cells <c>1cee505c</c> and <c>20730103</c>).
/// </summary>
/// <remarks>
/// <para>
/// Algorithm:
/// </para>
/// <list type="number">
/// <item>Build a peptide-to-proteins inverted index.</item>
/// <item>Subsumable proteins: any protein A whose peptides are a strict
/// subset of some other protein B's. The fast finder intersects
/// <c>peptideToProteins[pep]</c> over A's peptides; the intersection is
/// exactly the proteins that contain every one of A's peptides. The
/// smallest such B (lexicographic) is A's subsumer. O(N * avg|peps_A|) on
/// average, far cheaper than the O(N^2) PRISM reference on real datasets.</item>
/// <item>Collapse indistinguishable proteins (identical peptide sets) by
/// fingerprinting on the sorted-peptide string.</item>
/// <item>A peptide is "unique" to a canonical group iff exactly one
/// canonical group contains it after collapse; otherwise it is "razor" and
/// is greedy-assigned to the canonical group with the most unique peptides
/// (ties broken by lower canonical id). The notebook proves the iterative
/// PRISM reference is equivalent to a single sorted pass here because the
/// priority key is static across iterations.</item>
/// </list>
/// </remarks>
public static class ParsimonyEngine
{
    /// <summary>
    /// Maps each peptide to <c>(canonicalProteinGroup, "unique" | "razor")</c>.
    /// Peptides that map to zero active proteins after subsumption are absent
    /// from the result.
    /// </summary>
    public static Dictionary<string, (string Group, string Type)> Assign(
        IReadOnlyDictionary<string, IReadOnlyList<string>> peptideToProteins)
        => AssignWithAccessionGroups(peptideToProteins).Peptides;

    /// <summary>
    /// Same assignment as <see cref="Assign"/>, but also returns a
    /// protein-accession -> canonical-group map covering every protein in the
    /// input (active, indistinguishable, and subsumed). Used by the DIA-NN +
    /// Carafe backfill path to route predicted peptides into the SAME canonical
    /// group their confirmed protein landed in - without re-running parsimony
    /// over the combined evidence, which would perturb the confirmed grouping.
    /// </summary>
    public static (Dictionary<string, (string Group, string Type)> Peptides,
                   IReadOnlyDictionary<string, string> AccessionToGroup) AssignWithAccessionGroups(
        IReadOnlyDictionary<string, IReadOnlyList<string>> peptideToProteins)
    {
        // peptide -> set of proteins, protein -> set of peptides
        var pepToProts = new Dictionary<string, HashSet<string>>(peptideToProteins.Count);
        var protToPeps = new Dictionary<string, HashSet<string>>();
        foreach (var (peptide, prots) in peptideToProteins)
        {
            var set = new HashSet<string>(prots);
            pepToProts[peptide] = set;
            foreach (var prot in set)
            {
                if (!protToPeps.TryGetValue(prot, out var peps))
                {
                    peps = new HashSet<string>();
                    protToPeps[prot] = peps;
                }
                peps.Add(peptide);
            }
        }

        // Subsumable proteins via the inverted index.
        var subsumed = FindSubsumable(protToPeps);

        var active = new HashSet<string>(protToPeps.Keys.Where(p => !subsumed.ContainsKey(p)));

        // Collapse indistinguishable groups (identical peptide sets) by
        // fingerprinting on sorted-peptide string.
        var fingerprintToProteins = new Dictionary<string, List<string>>();
        foreach (var prot in active)
        {
            var key = SortedPeptideKey(protToPeps[prot]);
            if (!fingerprintToProteins.TryGetValue(key, out var list))
            {
                list = new List<string>();
                fingerprintToProteins[key] = list;
            }
            list.Add(prot);
        }
        var canonicalToPeps = new Dictionary<string, HashSet<string>>();
        foreach (var members in fingerprintToProteins.Values)
        {
            string canon = members.Min(StringComparer.Ordinal)!;
            canonicalToPeps[canon] = protToPeps[canon];
        }

        // peptide -> set of canonical groups that still claim it.
        var pepToCanons = new Dictionary<string, HashSet<string>>();
        foreach (var (canon, peps) in canonicalToPeps)
        {
            foreach (var pep in peps)
            {
                if (!pepToCanons.TryGetValue(pep, out var s))
                {
                    s = new HashSet<string>();
                    pepToCanons[pep] = s;
                }
                s.Add(canon);
            }
        }

        var uniquePerCanon = new Dictionary<string, HashSet<string>>();
        var sharedPeps = new HashSet<string>();
        foreach (var (pep, canons) in pepToCanons)
        {
            if (canons.Count == 1)
            {
                string canon = canons.First();
                if (!uniquePerCanon.TryGetValue(canon, out var set))
                {
                    set = new HashSet<string>();
                    uniquePerCanon[canon] = set;
                }
                set.Add(pep);
            }
            else
            {
                sharedPeps.Add(pep);
            }
        }

        // Razor assignment in static priority order: (most unique peptides,
        // then lowest canonical id). The PRISM iterative reference rescans
        // candidates each round; because the priority key is unchanged
        // across rounds a single sorted pass yields the same assignment.
        var priority = canonicalToPeps.Keys.OrderBy(c =>
        {
            int n = uniquePerCanon.TryGetValue(c, out var s) ? s.Count : 0;
            return -n;
        }).ThenBy(c => c, StringComparer.Ordinal).ToList();

        var razorPerCanon = new Dictionary<string, HashSet<string>>();
        var remaining = new HashSet<string>(sharedPeps);
        foreach (var canon in priority)
        {
            if (remaining.Count == 0) break;
            var canonPeps = canonicalToPeps[canon];
            var claim = new HashSet<string>();
            foreach (var pep in canonPeps)
            {
                if (remaining.Contains(pep)) claim.Add(pep);
            }
            if (claim.Count > 0)
            {
                razorPerCanon[canon] = claim;
                foreach (var pep in claim) remaining.Remove(pep);
            }
        }

        var result = new Dictionary<string, (string Group, string Type)>();
        foreach (var (canon, peps) in uniquePerCanon)
        {
            foreach (var pep in peps)
                result[pep] = (canon, "unique");
        }
        foreach (var (canon, peps) in razorPerCanon)
        {
            foreach (var pep in peps)
                result[pep] = (canon, "razor");
        }

        // Protein-accession -> canonical group, for EVERY input protein. Active
        // proteins map to the canon of their fingerprint (identical-peptide-set)
        // collapse; subsumed proteins follow the subsumption chain to an active
        // protein and inherit its canon. Same canon strings the peptide result
        // uses, so a predicted peptide routed by accession lands in exactly the
        // group its confirmed protein occupies.
        var accessionToGroup = new Dictionary<string, string>();
        foreach (var members in fingerprintToProteins.Values)
        {
            string canon = members.Min(StringComparer.Ordinal)!;
            foreach (var m in members) accessionToGroup[m] = canon;
        }
        foreach (var protA in subsumed.Keys)
        {
            string b = subsumed[protA];
            var guard = new HashSet<string> { protA };
            while (subsumed.TryGetValue(b, out var next) && guard.Add(b)) b = next;
            if (accessionToGroup.TryGetValue(b, out var canon))
                accessionToGroup[protA] = canon;
        }

        return (result, accessionToGroup);
    }

    /// <summary>
    /// For each protein A, find the smallest-id protein B whose peptide set
    /// is a strict superset of A's (i.e. B contains every peptide A has and
    /// at least one more). Returns <c>{ A: B }</c>. A is "subsumed by" B.
    /// </summary>
    internal static Dictionary<string, string> FindSubsumable(
        IReadOnlyDictionary<string, HashSet<string>> protToPeps)
    {
        // Build peptide -> set of containing proteins.
        var pepToProts = new Dictionary<string, HashSet<string>>();
        foreach (var (prot, peps) in protToPeps)
        {
            foreach (var pep in peps)
            {
                if (!pepToProts.TryGetValue(pep, out var s))
                {
                    s = new HashSet<string>();
                    pepToProts[pep] = s;
                }
                s.Add(prot);
            }
        }

        var subsumed = new Dictionary<string, string>();
        var proteinsSorted = protToPeps.Keys.OrderBy(p => p, StringComparer.Ordinal).ToList();
        foreach (var protA in proteinsSorted)
        {
            var pepsA = protToPeps[protA];
            if (pepsA.Count == 0) continue;

            // Pivot on the peptide with the smallest containing-protein set;
            // that minimises the candidate set we'll intersect with.
            string pivot = pepsA.First();
            int pivotSize = pepToProts[pivot].Count;
            foreach (var pep in pepsA)
            {
                int sz = pepToProts[pep].Count;
                if (sz < pivotSize)
                {
                    pivot = pep;
                    pivotSize = sz;
                }
            }
            var candidates = new HashSet<string>(pepToProts[pivot]);
            candidates.Remove(protA);

            foreach (var pep in pepsA)
            {
                if (pep == pivot) continue;
                candidates.IntersectWith(pepToProts[pep]);
                if (candidates.Count == 0) break;
            }
            if (candidates.Count == 0) continue;

            // Smallest-id candidate with a strict superset of A's peptides
            // (i.e. strictly larger peptide count) is A's subsumer. Same-size
            // candidates are "indistinguishable" - handled by the fingerprint
            // collapse step in Assign.
            foreach (var protB in candidates.OrderBy(c => c, StringComparer.Ordinal))
            {
                if (protToPeps[protB].Count > pepsA.Count)
                {
                    subsumed[protA] = protB;
                    break;
                }
            }
        }
        return subsumed;
    }

    private static string SortedPeptideKey(HashSet<string> peps)
    {
        var sorted = peps.ToArray();
        Array.Sort(sorted, StringComparer.Ordinal);
        return string.Join("", sorted);
    }
}
