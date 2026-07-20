using System.Globalization;
using System.Text;
using SkylineCadenza.Core.Scheduling;

namespace SkylineCadenza.Core.Output;

/// <summary>
/// Writes a peptide-style transition list CSV in the column layout
/// Skyline's <c>--import-transition-list=</c> expects for proteomics
/// documents. Used by the push-to-Skyline flow to populate the target
/// tree after the assay BLIB has been written via
/// <see cref="BlibAssayWriter"/> and registered as a peptide-settings
/// library; the BLIB carries spectra and per-replicate RT boundaries
/// for chromatogram extraction, and this transition list tells Skyline
/// exactly which (peptide, charge, fragment) triples to display in the
/// document tree.
/// </summary>
/// <remarks>
/// The builder emits one row per (scheduled precursor, top-N library
/// fragment). Each fragment row carries its real <c>Product Charge</c>
/// from <see cref="FragmentIon.Charge"/> - earlier versions hardcoded
/// <c>Product Charge = 1</c>, which produced "no matching product ion"
/// errors when DIA-NN-picked fragments were actually +2 (common for
/// longer tryptic peptides).
///
/// Modification syntax is normalised to Skyline's preferred
/// <c>C[UniMod:4]</c> / <c>C[+57.0]</c> form: brackets only, no
/// flanking underscores, no parentheses.
/// </remarks>
public static class PeptideTransitionListBuilder
{
    /// <summary>Number of fragment rows emitted per peptide.</summary>
    public const int FragmentsPerPeptide = BlibAssayWriter.PeaksPerSpectrum;

    public static string Build(
        IReadOnlyList<Candidate> candidates,
        ScheduleResult schedule)
    {
        var sb = new StringBuilder();
        sb.AppendLine(string.Join(",",
            "Protein Name",
            "Peptide Modified Sequence",
            "Precursor Charge",
            "Precursor m/z",
            "Product m/z",
            "Product Charge",
            "Explicit Retention Time",
            "Explicit Retention Time Window",
            "Note"));

        var inv = CultureInfo.InvariantCulture;
        var slotsById = schedule.Slots.ToDictionary(s => s.Id);

        for (int i = 0; i < schedule.ScheduledIndices.Length; i++)
        {
            int candIdx = schedule.ScheduledIndices[i];
            int slotId = schedule.ScheduledSlotIds[i];
            var c = candidates[candIdx];

            // Skip heavy isotope-labeled precursors (e.g. spiked PRTC
            // standards). Skyline's headless transition-list importer
            // rejects an explicit isotope label it can't match to a
            // registered modification ("isotope labeled entry found
            // without matching settings"), and because such a row is
            // used to validate the precursor-m/z column, a single one
            // aborts the whole import with the misleading "No valid
            // precursor m/z column found, line 1". These precursors are
            // already added to the document as heavy through the
            // registered spectral library, so we omit them here.
            if (!c.IsLightLabel) continue;

            var slot = slotsById[slotId];

            double windowMin = slot.RtStop - slot.RtStart;
            string note = $"slot={slotId};type={c.PeptideType}";
            string protein = c.ProteinGroup;
            string peptide = NormalizeModifiedSequence(c.ModifiedSequence);

            // Precursor isotope rows first (M+0, M+1, M+2) so they
            // appear as precursor transitions in the document tree
            // alongside the fragment rows. Skyline recognizes
            // Product Charge == Precursor Charge with Product m/z
            // close to the precursor's monoisotope spacing as a
            // precursor transition. Mass step = 1.003355 Da per ^13C
            // -> ^12C neutron offset.
            for (int isotope = 0; isotope < PrecursorIsotopesPerPeptide; isotope++)
            {
                double isotopeMz = c.PrecursorMz
                    + NeutronMassDa * isotope / c.PrecursorCharge;
                AppendRow(sb, inv, protein, peptide, c.PrecursorCharge,
                    c.PrecursorMz,
                    productMz: isotopeMz,
                    productCharge: c.PrecursorCharge,
                    c.RtApex, windowMin, note);
            }

            if (c.Fragments.Length == 0)
            {
                AppendRow(sb, inv, protein, peptide, c.PrecursorCharge,
                    c.PrecursorMz, productMz: null, productCharge: null,
                    c.RtApex, windowMin, note);
                continue;
            }

            int take = Math.Min(FragmentsPerPeptide, c.Fragments.Length);
            for (int k = 0; k < take; k++)
            {
                var frag = c.Fragments[k];
                AppendRow(sb, inv, protein, peptide, c.PrecursorCharge,
                    c.PrecursorMz,
                    productMz: frag.Mz,
                    productCharge: frag.Charge,
                    c.RtApex, windowMin, note);
            }
        }
        return sb.ToString();
    }

    /// <summary>
    /// Mass difference between two adjacent precursor isotopes in
    /// Daltons. The monoisotopic envelope's spacing in m/z space is
    /// this value divided by the precursor charge.
    /// </summary>
    private const double NeutronMassDa = 1.003355;

    /// <summary>
    /// Number of precursor isotopes (M+0, M+1, M+2) emitted per
    /// peptide. Must match the
    /// <c>--full-scan-precursor-threshold</c> value
    /// <see cref="SkylineRpc.SkylineSettingsConfigurator"/> sends so
    /// the document tree's precursor nodes line up with what
    /// Skyline's MS1 chromatogram extractor actually pulls.
    /// </summary>
    public const int PrecursorIsotopesPerPeptide = 3;

    private static void AppendRow(
        StringBuilder sb, CultureInfo inv,
        string protein, string peptide, int precursorCharge,
        double precursorMz, double? productMz, int? productCharge,
        double rtApex, double windowMin, string note)
    {
        sb.AppendLine(string.Join(",",
            Csv(protein),
            Csv(peptide),
            precursorCharge.ToString(inv),
            precursorMz.ToString("0.0000", inv),
            productMz.HasValue ? productMz.Value.ToString("0.0000", inv) : "",
            productCharge.HasValue ? productCharge.Value.ToString(inv) : "",
            rtApex.ToString("0.0000", inv),
            windowMin.ToString("0.0000", inv),
            Csv(note)));
    }

    /// <summary>
    /// Convert a DIA-NN <c>C(UniMod:4)</c> or Carafe <c>_C[UniMod:4]_</c>
    /// style modified sequence to Skyline's preferred <c>C[UniMod:4]</c>
    /// form. Strips flanking underscores and swaps parentheses for square
    /// brackets.
    /// </summary>
    private static string NormalizeModifiedSequence(string modSeq)
    {
        if (string.IsNullOrEmpty(modSeq)) return modSeq;
        var sb = new StringBuilder(modSeq.Length);
        foreach (var ch in modSeq)
        {
            switch (ch)
            {
                case '_': break;
                case '(': sb.Append('['); break;
                case ')': sb.Append(']'); break;
                default: sb.Append(ch); break;
            }
        }
        return sb.ToString();
    }

    private static string Csv(string value)
    {
        if (value.IndexOfAny(new[] { ',', '"', '\n' }) < 0) return value;
        return "\"" + value.Replace("\"", "\"\"") + "\"";
    }
}
