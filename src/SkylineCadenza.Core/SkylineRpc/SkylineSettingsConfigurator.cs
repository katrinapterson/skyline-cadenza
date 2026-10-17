#nullable enable

using System.Text;
using SkylineCadenza.Core.Output;
using SkylineCadenza.Core.Scheduling;

namespace SkylineCadenza.Core.SkylineRpc;

/// <summary>
/// Aligns the running Skyline document's peptide / transition-filter
/// settings with the shape of the assay Cadenza is about to push.
/// </summary>
/// <remarks>
/// <para>
/// All updates go through <see cref="IJsonToolService.RunCommand"/>.
/// Per Nick Shulman: anything <c>SkylineCmd</c> can do is also
/// reachable via JSON-RPC <c>RunCommand</c> against the live document,
/// so we use the same flags here that the command-line build accepts.
/// Some flags may not be recognised by every Skyline version - those
/// surface as text in the returned report rather than aborting the
/// rest of the configuration.
/// </para>
/// <para>
/// The recommendation is intentionally union-of-assay for charge sets
/// and the broadest sensible ion-type / ion-range filter, so the
/// user's downstream Skyline analysis isn't surprised by a narrower
/// filter than what their library actually contains.
/// </para>
/// </remarks>
public static class SkylineSettingsConfigurator
{
    /// <summary>
    /// Recommended document settings derived from the scheduled assay.
    /// </summary>
    public sealed record Recommendation(
        IReadOnlyList<int> PrecursorIonCharges,
        IReadOnlyList<int> ProductIonCharges,
        IReadOnlyList<string> ProductIonTypes,
        int LibraryPickTopN,
        int PeptideMinLength,
        int PeptideMaxLength,
        int PeptideExcludeNTerminalAAs,
        string FullScanAcquisitionMethod,
        string? FullScanIsolationScheme,
        string RetentionTimeFilter,
        double RetentionTimeFilterToleranceMin,
        string PrecursorIsotopes,
        int PrecursorIsotopeCount,
        string PrecursorAnalyzer,
        double PrecursorMassAccuracyPpm)
    {
        /// <summary>
        /// Short single-line description for the status bar.
        /// </summary>
        public string ToStatusLine() =>
            $"acquisition method {FullScanAcquisitionMethod}, "
            + (FullScanIsolationScheme is null
                ? ""
                : $"isolation scheme '{FullScanIsolationScheme}', ")
            + $"MS1 isotopes {PrecursorIsotopeCount} (M+0..M+{PrecursorIsotopeCount - 1}) "
            + $"@ {PrecursorAnalyzer} {PrecursorMassAccuracyPpm:0} ppm, "
            + $"RT filter {RetentionTimeFilter} ±{RetentionTimeFilterToleranceMin:0.0} min, "
            + $"precursor charges {{{string.Join(",", PrecursorIonCharges)}}}, "
            + $"product charges {{{string.Join(",", ProductIonCharges)}}}, "
            + $"ion types {{{string.Join(",", ProductIonTypes)}}}, "
            + $"library pick top {LibraryPickTopN}, "
            + $"peptide length {PeptideMinLength}-{PeptideMaxLength}, "
            + $"exclude {PeptideExcludeNTerminalAAs} N-terminal AAs.";

        /// <summary>
        /// Step-by-step instructions the user can follow in the Skyline UI
        /// to align the document with the recommendation (used as a
        /// fallback when one of the <c>SkylineCmd</c> flags fails).
        /// </summary>
        public string ToUiInstructions()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Verify in Skyline (Settings menu):");
            sb.AppendLine("  Transition Settings > Full-Scan:");
            sb.AppendLine($"    Acquisition method: {FullScanAcquisitionMethod}");
            sb.AppendLine("  Transition Settings > Filter > Peptides:");
            sb.AppendLine($"    Precursor charges: {string.Join(", ", PrecursorIonCharges)}");
            sb.AppendLine($"    Ion charges:      {string.Join(", ", ProductIonCharges)}");
            sb.AppendLine($"    Ion types:        {string.Join(", ", ProductIonTypes)}");
            sb.AppendLine($"    Product ion selection: ion 1 to last ion");
            sb.AppendLine("  Transition Settings > Library:");
            sb.AppendLine($"    Pick {LibraryPickTopN} product ions");
            sb.AppendLine("    Match tolerance m/z: 0.5 (or your standard)");
            sb.AppendLine("  Peptide Settings > Filter:");
            sb.AppendLine($"    Min length: {PeptideMinLength}");
            sb.AppendLine($"    Max length: {PeptideMaxLength}");
            sb.AppendLine($"    Exclude N-terminal AAs: {PeptideExcludeNTerminalAAs}");
            return sb.ToString().TrimEnd();
        }
    }

    /// <summary>
    /// Compute the recommendation from the scheduled candidates, the
    /// acquisition mode, and the firing-pad value used at scheduling
    /// time.
    /// <list type="bullet">
    /// <item>PRM mode -&gt; Skyline's <c>PRM</c> acquisition method, no
    /// isolation scheme (each precursor has its own quadrupole window).</item>
    /// <item>MTM mode -&gt; <c>DIA</c> acquisition method + isolation
    /// scheme <c>"Results only"</c>. MTM is multiplexed-targeted DIA
    /// from Skyline's perspective (each MS/MS spectrum carries multiple
    /// co-isolated precursors), and we don't pre-specify isolation
    /// windows because Cadenza schedules them per slot - Skyline reads
    /// the actual isolation scheme from the raw file at import time.</item>
    /// </list>
    /// The RT filter is <c>scheduling_windows</c> ("Use only scans
    /// within X minutes of predicted RT"). The "predicted RT" comes
    /// from the document's <c>PeptideSettings.Prediction</c>, which
    /// the user must configure manually one time (no SkylineCmd flag
    /// exposes the predictor or the "Use measured retention times
    /// when present" toggle): set the alignment target to library RTs
    /// so the predicted apex falls back to the BLIB's library RT for
    /// each peptide. We use this in preference to <c>ms2_ids</c>
    /// because freshly-acquired PRM / MTM raw files don't carry any
    /// MS/MS IDs - the assay was DESIGNED by Cadenza, not searched -
    /// so <c>ms2_ids</c> would have nothing to filter against once
    /// the data are imported. The tolerance is the assay's worst-case
    /// half-width from apex to peak edge, padded by
    /// <paramref name="firingPadMin"/>, rounded up to 0.1 min.
    /// </summary>
    public static Recommendation Recommend(
        IReadOnlyList<Candidate> scheduledCandidates,
        AcquisitionMode mode,
        double firingPadMin = 0.25)
    {
        var precursorCharges = scheduledCandidates
            .Select(c => c.PrecursorCharge)
            .Where(z => z > 0)
            .Distinct()
            .OrderBy(z => z)
            .ToList();
        if (precursorCharges.Count == 0) precursorCharges.AddRange(new[] { 2, 3 });

        var fragmentCharges = scheduledCandidates
            .SelectMany(c => c.Fragments.Select(f => f.Charge))
            .Where(z => z > 0)
            .Distinct()
            .OrderBy(z => z)
            .ToList();
        if (fragmentCharges.Count == 0) fragmentCharges.Add(1);

        var lengths = scheduledCandidates
            .Select(c => c.StrippedSequence?.Length ?? 0)
            .Where(n => n > 0)
            .ToList();
        int minLen = lengths.Count == 0 ? 6 : Math.Max(4, lengths.Min());
        int maxLen = lengths.Count == 0 ? 30 : Math.Min(60, lengths.Max());

        string acquisitionMethod = mode == AcquisitionMode.Prm ? "PRM" : "DIA";
        string? isolationScheme = mode == AcquisitionMode.Prm ? null : "Results only";

        // RT filter tolerance: worst-case half-width from apex to peak
        // edge across the scheduled assay, plus the firing pad, rounded
        // up to 0.1 min. The floor of 0.2 min keeps the tolerance large
        // enough to be useful even when every peptide has a near-zero
        // measured peak width (typically a synthesized fallback).
        double maxHalfWidth = 0;
        foreach (var c in scheduledCandidates)
        {
            double belowApex = c.RtApex - c.RtStart;
            double aboveApex = c.RtStop - c.RtApex;
            double half = Math.Max(belowApex, aboveApex);
            if (half > maxHalfWidth) maxHalfWidth = half;
        }
        double rawTolerance = maxHalfWidth + firingPadMin;
        double rtFilterToleranceMin = Math.Max(0.2, Math.Ceiling(rawTolerance * 10) / 10);

        return new Recommendation(
            PrecursorIonCharges: precursorCharges,
            ProductIonCharges: fragmentCharges,
            // y/b are the assay's actual fragment families; we also
            // include p (precursor) so Skyline keeps the precursor m/z
            // tracked as an MS1 reference transition. The transition
            // list Cadenza imports only carries y/b rows, so adding p
            // here is additive for the document's transition picker -
            // it doesn't change anything we push.
            ProductIonTypes: new[] { "y", "b", "p" },
            LibraryPickTopN: BlibAssayWriter.PeaksPerSpectrum,
            PeptideMinLength: minLen,
            PeptideMaxLength: maxLen,
            FullScanAcquisitionMethod: acquisitionMethod,
            FullScanIsolationScheme: isolationScheme,
            RetentionTimeFilter: "scheduling_windows",
            RetentionTimeFilterToleranceMin: rtFilterToleranceMin,
            // MS1 filtering: top 3 precursor isotopes (M+0, M+1, M+2)
            // on a centroided analyzer at 10 ppm. Matches the standard
            // PRM / MTM workflow on Orbitrap / Stellar; these values
            // came from Mike's working Transition Settings > Full-Scan
            // dialog.
            PrecursorIsotopes: "Count",
            PrecursorIsotopeCount: 3,
            PrecursorAnalyzer: "centroided",
            PrecursorMassAccuracyPpm: 10,
            // 0 = include every position. Skyline's default of 25
            // (intended to skip signal-peptide N-terminal fragments)
            // would drop peptides Cadenza explicitly designed the
            // assay around, since we picked them from observation
            // not from a sequence-position filter.
            PeptideExcludeNTerminalAAs: 0);
    }

    /// <summary>
    /// Result of running the configurator: the recommendation that was
    /// computed, plus any <c>SkylineCmd</c> output that came back from
    /// the live document (empty / whitespace means Skyline accepted the
    /// settings silently).
    /// </summary>
    public sealed record ConfigureResult(Recommendation Recommendation, string SkylineOutput);

    /// <summary>
    /// The SkylineCmd flags that apply <paramref name="rec"/>, in the
    /// batches <see cref="ConfigureAsync"/> sends them in. Shared with
    /// <see cref="Output.SkylineDocumentExport"/>, so a document built
    /// headless gets the same settings as a push to a running Skyline.
    /// </summary>
    public static string[][] ArgumentBatches(Recommendation rec)
    {
        // One RunCommand per individual flag. Verified against
        // pwiz_tools/Skyline/CommandArgs.cs:
        //
        //   --pep-min-length / --pep-max-length          (line 1616, 1618)
        //   --tran-precursor-ion-charges                 (... look near 1530s)
        //   --tran-product-ion-charges                   (line 1537)
        //   --tran-product-ion-types                     (line 1540)
        //   --tran-product-start-ion                     (line 1542)
        //   --tran-product-end-ion                       (line 1545)
        //   --library-pick-product-ions=filter           (line 1564, enum: filter)
        //   --library-product-ions=N                     (line 1560)
        //
        // Batches are kept narrow because Skyline's SkylineCmd parser
        // aborts the entire batch on the first unknown argument (it
        // prints "Error: Unexpected argument X" and exits). One flag
        // per batch means any future unknown flag we add by mistake
        // only loses its own setting, not the ones around it.
        var batchList = new List<string[]>
        {
            new[] { "--pep-min-length=" + rec.PeptideMinLength },
            new[] { "--pep-max-length=" + rec.PeptideMaxLength },
            new[] { "--pep-exclude-nterminal-aas=" + rec.PeptideExcludeNTerminalAAs },
            new[] { "--tran-precursor-ion-charges=" + string.Join(",", rec.PrecursorIonCharges) },
            new[] { "--tran-product-ion-charges=" + string.Join(",", rec.ProductIonCharges) },
            new[] { "--tran-product-ion-types=" + string.Join(",", rec.ProductIonTypes) },
            new[] { "--tran-product-start-ion=ion 1" },
            new[] { "--tran-product-end-ion=last ion" },
            new[] { "--library-pick-product-ions=filter" },
            new[] { "--library-product-ions=" + rec.LibraryPickTopN },
        };
        // DIA mode + isolation scheme are mutually constrained and must
        // be sent in a single RunCommand call so Skyline validates the
        // end state, not each intermediate. Splitting them deadlocks:
        // - Setting --full-scan-isolation-scheme alone while the
        //   document is still in AcquisitionMethod = None throws
        //   "No other full-scan MS/MS filter settings are allowed
        //   when precursor filter is none" (rejected silently with no
        //   stdout echo).
        // - Setting --full-scan-acquisition-method=DIA alone while
        //   IsolationScheme is still null throws "An isolation window
        //   width value is required in DIA mode."
        // Bundling them lets Skyline apply both fields and run
        // DoValidate once against the final (DIA + Results only)
        // state which is the valid configuration. PRM mode has no
        // isolation scheme so we just send the acquisition method.
        if (rec.FullScanIsolationScheme is not null)
        {
            batchList.Add(new[]
            {
                "--full-scan-isolation-scheme=" + rec.FullScanIsolationScheme,
                "--full-scan-acquisition-method=" + rec.FullScanAcquisitionMethod,
            });
        }
        else
        {
            batchList.Add(new[] { "--full-scan-acquisition-method=" + rec.FullScanAcquisitionMethod });
        }
        // RT filter: tell Skyline to extract chromatograms within
        // the assay's worst-case half-width of each peptide's
        // BLIB-derived predicted apex.
        batchList.Add(new[] { "--full-scan-rt-filter=" + rec.RetentionTimeFilter });
        batchList.Add(new[] { "--full-scan-rt-filter-tolerance="
            + rec.RetentionTimeFilterToleranceMin.ToString(
                "0.0", System.Globalization.CultureInfo.InvariantCulture) });
        // MS1 isotope filter: top N precursor isotopes on a centroided
        // analyzer at 10 ppm. Bundled into a single RunCommand call
        // because Skyline's TransitionFullScan.DoValidate runs against
        // the combined end state: PrecursorIsotopes != None requires
        // the analyzer + resolution to be configured at the same time
        // (ValidateRes(PrecursorMassAnalyzer, PrecursorRes,
        // PrecursorResMz) at TransitionSettings.cs line 2881). Setting
        // them individually trips the same kind of mutual-lockout we
        // hit on DIA + Isolation scheme.
        batchList.Add(new[]
        {
            "--full-scan-precursor-isotopes=" + rec.PrecursorIsotopes,
            "--full-scan-precursor-threshold=" + rec.PrecursorIsotopeCount.ToString(
                System.Globalization.CultureInfo.InvariantCulture),
            "--full-scan-precursor-analyzer=" + rec.PrecursorAnalyzer,
            "--full-scan-precursor-res=" + rec.PrecursorMassAccuracyPpm.ToString(
                "0.#", System.Globalization.CultureInfo.InvariantCulture),
        });
        return batchList.ToArray();
    }

    /// <summary>
    /// Compute the recommendation and apply it to the live document via
    /// <c>RunCommand</c>. Per Nick Shulman, <c>RunCommand</c> against
    /// the JSON-RPC pipe operates on the document held by the running
    /// Skyline UI, so the same flags <c>SkylineCmd</c> would accept on
    /// a saved file also work here.
    /// </summary>
    public static async Task<ConfigureResult> ConfigureAsync(
        SkylineSession session,
        IReadOnlyList<Candidate> scheduledCandidates,
        AcquisitionMode mode,
        double firingPadMin = 0.25,
        CancellationToken cancellationToken = default)
    {
        var rec = Recommend(scheduledCandidates, mode, firingPadMin);

        string[][] argBatches = ArgumentBatches(rec);

        var collected = new List<string>(argBatches.Length);
        foreach (var args in argBatches)
        {
            cancellationToken.ThrowIfCancellationRequested();
            try
            {
                string output = await Task.Run(
                    () => session.Execute(c => c.RunCommand(args)),
                    cancellationToken);
                if (!string.IsNullOrWhiteSpace(output))
                    collected.Add(output.Trim());
            }
            catch (Exception ex)
            {
                collected.Add($"({string.Join(" ", args)}): {ex.Message}");
            }
        }
        return new ConfigureResult(rec, string.Join(" | ", collected));
    }
}
