namespace SkylineCadenza.Core.Scheduling;

/// <summary>
/// One scheduling candidate. Represents a single precursor (peptide + charge)
/// already assigned to a parsimonious protein group with its top-4 fragment
/// m/z values from the spectral library.
/// </summary>
/// <remarks>
/// Mirrors the row layout of the <c>candidates</c> DataFrame in the
/// targeted-modeling notebook (cell <c>20730103</c>). The unpadded
/// <see cref="RtStart"/> / <see cref="RtStop"/> drive the co-elution
/// intersection check; the scheduler pads them by
/// <see cref="SchedulingParameters.FiringPadMin"/> when computing the slot's
/// firing window for the budget.
/// </remarks>
public sealed class Candidate
{
    public required string PrecursorId { get; init; }
    public required string StrippedSequence { get; init; }
    public required string ModifiedSequence { get; init; }
    public required int PrecursorCharge { get; init; }
    public required double PrecursorMz { get; init; }

    /// <summary>Unpadded peak start (min).</summary>
    public required double RtStart { get; init; }

    /// <summary>Unpadded peak stop (min).</summary>
    public required double RtStop { get; init; }

    /// <summary>Apex retention time (min) - written into the Skyline target list.</summary>
    public required double RtApex { get; init; }

    public required double PrecursorQuantity { get; init; }

    /// <summary>
    /// Predicted within-protein priority (higher = better). Populated from a
    /// predicted-ranking overlay (Bromo <c>mean_pred_score</c> or Carafe
    /// <c>detectability_pred</c>) in Carafe-only mode. Null when no overlay was
    /// supplied. Consumed only by <see cref="PeptidePriority.PredictedRanking"/>
    /// for ordering peptides within a protein; it is NOT a cross-protein signal.
    /// </summary>
    public double? PredictedPriority { get; init; }

    /// <summary>Peptide-precursor q-value (DIA-NN <c>Q.Value</c>).</summary>
    public required double QValue { get; init; }

    /// <summary>Protein-group q-value (DIA-NN <c>PG.Q.Value</c>) or NaN if unavailable.</summary>
    public required double ProteinQValue { get; init; }

    /// <summary>1 if the peptide is unique to a single protein in the database.</summary>
    public required int Proteotypic { get; init; }

    /// <summary>Canonical protein-group id after parsimony (e.g. UniProt accession).</summary>
    public required string ProteinGroup { get; init; }

    /// <summary><c>"unique"</c> or <c>"razor"</c> per the parsimony assignment.</summary>
    public required string PeptideType { get; init; }

    /// <summary>
    /// Full predicted-fragment record kept on the candidate: up to
    /// <c>FragmentLimit</c> entries (intensity descending), with the m/z
    /// values, library intensities, and per-fragment charges from the
    /// originating ingest path. Used to write a BLIB at assay-push time
    /// and to derive <see cref="Top4Fragments"/> for the scheduler's
    /// fragment-clash check. May be empty for precursors with no library
    /// entry.
    /// </summary>
    public required FragmentIon[] Fragments { get; init; }

    /// <summary>
    /// Top-4 predicted fragment m/z values (sorted ascending) used by
    /// the scheduler's fragment-clash check. Set by the ingest path as
    /// the top-4 entries from <see cref="Fragments"/> by intensity,
    /// projected to m/z and re-sorted by m/z. Cached because the
    /// scheduler reads this in its inner loop.
    /// </summary>
    public required double[] Top4Fragments { get; init; }

    /// <summary>DIA-NN <c>Run</c> name. Used to split per-run / per-GPF views.</summary>
    public required string Run { get; init; }

    /// <summary>
    /// Number of fragments Cadenza keeps per precursor. Set wide enough
    /// to back both the top-4 clash check and the top-6 written into the
    /// assay BLIB without needing to re-ingest when the BLIB writer's
    /// width changes.
    /// </summary>
    public const int FragmentLimit = 12;

    /// <summary>
    /// Project the top-k entries of <paramref name="fragments"/> (already
    /// intensity-sorted descending) onto a sorted-ascending m/z array,
    /// matching <see cref="Top4Fragments"/>'s contract.
    /// </summary>
    public static double[] DeriveTopMz(IReadOnlyList<FragmentIon> fragments, int k = 4)
    {
        int take = Math.Min(k, fragments.Count);
        if (take == 0) return Array.Empty<double>();
        var mz = new double[take];
        for (int i = 0; i < take; i++) mz[i] = fragments[i].Mz;
        Array.Sort(mz);
        return mz;
    }
}
