#nullable enable

using System.Diagnostics;
using System.Text;
using SkylineCadenza.Core.SkylineRpc;

namespace SkylineCadenza.Core.Output;

/// <summary>
/// Writes a scheduled assay as a Skyline document on disk, without a running
/// Skyline: the same three pieces "Push to Skyline" sends to the live
/// document (the assay BLIB as the document's library, the settings
/// <see cref="SkylineSettingsConfigurator"/> recommends, and the peptide
/// transition list), applied to a new document with <c>SkylineCmd</c>. Used
/// by <c>SkylineCadenzaCli export-thermo ... --skyline OUT.sky</c>, so the
/// document matches the Thermo CSV written from the same schedule and the
/// targeted runs can be imported into it.
/// </summary>
public static class SkylineDocumentExport
{
    /// <param name="Document">The .sky written.</param>
    /// <param name="Library">The assay BLIB beside it, registered as the document's library.</param>
    /// <param name="TransitionList">The transition list imported into it, kept beside it.</param>
    /// <param name="Precursors">Precursors written to the library.</param>
    /// <param name="SkylineOutput">What SkylineCmd printed, for the log.</param>
    public sealed record Result(string Document, string Library, string TransitionList, int Precursors, string SkylineOutput);

    /// <summary>
    /// Writes the assay in <paramref name="export"/> to <paramref name="documentPath"/>
    /// (replacing it), with its BLIB and transition list beside it.
    /// </summary>
    /// <param name="skylineCmd">SkylineCmd.exe; null for the newest installed.</param>
    public static async Task<Result> WriteAsync(
        DiannThermoExport.Result export,
        string documentPath,
        string? skylineCmd = null,
        CancellationToken cancellationToken = default)
    {
        var cmd = skylineCmd ?? FindSkylineCmd()
            ?? throw new InvalidOperationException("SkylineCmd was not found. Install Skyline, or pass --skylinecmd.");
        documentPath = Path.GetFullPath(documentPath);
        var folder = Path.GetDirectoryName(documentPath)!;
        Directory.CreateDirectory(folder);
        var name = Path.GetFileNameWithoutExtension(documentPath);

        // Drivers only, as the app's push does by default: passengers ride in
        // windows the instrument already fires, so they add extraction targets,
        // not acquisitions, and stay a choice made in the GUI.
        var blib = Path.Combine(folder, name + ".blib");
        var written = BlibAssayWriter.Write(export.Candidates, export.Schedule, blib, name);
        var transitionList = Path.Combine(folder, name + "_transitions.csv");
        await File.WriteAllTextAsync(transitionList,
            PeptideTransitionListBuilder.Build(export.Candidates, export.Schedule), cancellationToken);

        var scheduled = export.Schedule.ScheduledIndices.Select(i => export.Candidates[i]).ToList();
        var recommendation = SkylineSettingsConfigurator.Recommend(
            scheduled, export.Parameters.Mode, firingPadMin: export.Parameters.FiringPadSec / 60.0);

        var output = await RunAsync(cmd, Arguments(documentPath, blib, name, transitionList, recommendation), folder, cancellationToken);
        if (!File.Exists(documentPath))
            throw new InvalidOperationException($"SkylineCmd did not write {documentPath}: {Last(output)}");
        return new Result(documentPath, blib, transitionList, written.RefSpectraWritten, output);
    }

    /// <summary>
    /// One SkylineCmd call: a new document, the recommended settings, the
    /// library, then the transition list. Settings come before the import so
    /// the transition filter accepts the assay's +2 fragments.
    /// </summary>
    public static List<string> Arguments(
        string documentPath, string blib, string libraryName, string transitionList,
        SkylineSettingsConfigurator.Recommendation recommendation)
    {
        var args = new List<string> { "--new=" + documentPath, "--overwrite" };
        // The same flags as the push, including the "scheduling_windows" RT
        // filter: it needs no predictor here, because the transition list
        // gives every precursor an explicit retention time and window
        // (observed 2026-10-17 with Skyline-daily: SkylineCmd accepts it and
        // imports with it). Extracting the whole run instead would let MTM
        // windows at other times, which can hold the same m/z, offer peaks.
        foreach (var batch in SkylineSettingsConfigurator.ArgumentBatches(recommendation))
            args.AddRange(batch);
        args.Add("--add-library-path=" + blib);
        args.Add("--add-library-name=" + libraryName);
        args.Add("--import-transition-list=" + transitionList);
        // As the push: one rejected row must not discard the whole assay.
        args.Add("--ignore-transition-errors");
        args.Add("--save");
        return args;
    }

    /// <summary>
    /// The newest SkylineCmd.exe: Program Files\Skyline, or a ClickOnce
    /// install (Skyline or Skyline-daily) under %LOCALAPPDATA%\Apps\2.0,
    /// dated by the Skyline executable beside it.
    /// </summary>
    public static string? FindSkylineCmd()
    {
        var found = new List<(string Cmd, DateTime When)>();
        foreach (var root in new[] { Environment.SpecialFolder.ProgramFiles, Environment.SpecialFolder.ProgramFilesX86 })
        {
            var path = Path.Combine(Environment.GetFolderPath(root), "Skyline", "SkylineCmd.exe");
            if (File.Exists(path)) found.Add((path, File.GetLastWriteTimeUtc(path)));
        }

        var clickOnce = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "Apps", "2.0");
        if (Directory.Exists(clickOnce))
        {
            try
            {
                foreach (var cmd in Directory.EnumerateFiles(clickOnce, "SkylineCmd.exe",
                             new EnumerationOptions { RecurseSubdirectories = true, IgnoreInaccessible = true }))
                {
                    var folder = Path.GetDirectoryName(cmd)!;
                    var main = new[] { "Skyline-daily.exe", "Skyline.exe" }.Select(n => Path.Combine(folder, n)).FirstOrDefault(File.Exists);
                    if (main is not null) found.Add((cmd, File.GetLastWriteTimeUtc(main)));
                }
            }
            catch (IOException)
            {
                // Use what was found.
            }
        }

        return found.OrderByDescending(f => f.When).Select(f => f.Cmd).FirstOrDefault();
    }

    private static async Task<string> RunAsync(string cmd, IReadOnlyList<string> args, string workingDirectory, CancellationToken cancellationToken)
    {
        var info = new ProcessStartInfo(cmd)
        {
            UseShellExecute = false,
            CreateNoWindow = true,
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            StandardOutputEncoding = Encoding.UTF8,
            StandardErrorEncoding = Encoding.UTF8,
            WorkingDirectory = workingDirectory,
        };
        foreach (var arg in args) info.ArgumentList.Add(arg);

        using var process = Process.Start(info) ?? throw new InvalidOperationException("SkylineCmd could not be started.");
        var stdout = process.StandardOutput.ReadToEndAsync(CancellationToken.None);
        var stderr = process.StandardError.ReadToEndAsync(CancellationToken.None);
        try
        {
            await process.WaitForExitAsync(cancellationToken);
        }
        catch (OperationCanceledException)
        {
            try { process.Kill(entireProcessTree: true); } catch (InvalidOperationException) { }
            throw;
        }

        var text = await stdout + await stderr;
        var errors = text.Split('\n').Select(l => l.Trim())
            .Where(l => l.StartsWith("Error", StringComparison.OrdinalIgnoreCase)).ToList();
        if (process.ExitCode != 0 || errors.Count > 0)
            throw new InvalidOperationException(
                "SkylineCmd could not build the assay document: " + (errors.Count > 0 ? string.Join(" ", errors.Take(3)) : Last(text)));
        return text;
    }

    private static string Last(string text) =>
        text.Trim().Split('\n').Select(l => l.Trim()).LastOrDefault(l => l.Length > 0) ?? "(no output)";
}
