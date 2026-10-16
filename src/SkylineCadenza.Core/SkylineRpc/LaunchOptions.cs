namespace SkylineCadenza.Core.SkylineRpc;

/// <summary>
/// The app's command line. Two shapes:
/// <list type="bullet">
/// <item>Launched from Skyline's Tools menu: <c>args[0]</c> is the Skyline
/// pipe name (<c>$(SkylineConnection)</c>), or no arguments at all (attach to
/// the most recent running Skyline). Handled exactly as before.</item>
/// <item>Standalone, driven by another app: the first argument starts with
/// <c>--</c>. <c>--report &lt;report.parquet&gt;</c> preloads a DIA-NN report
/// and runs Load + schedule; <c>--export-dir &lt;folder&gt;</c> makes "Export
/// Thermo CSV..." default to that folder. No Skyline connection is made.</item>
/// </list>
/// </summary>
public sealed record LaunchOptions
{
    /// <summary>True when the first argument is an option, i.e. not a Skyline pipe name.</summary>
    public bool Standalone { get; init; }

    /// <summary>Arguments handed to <see cref="SkylineSession.FromArguments"/> when not standalone.</summary>
    public string[] SkylineArgs { get; init; } = Array.Empty<string>();

    /// <summary>DIA-NN report to preload (<c>--report</c>), or null.</summary>
    public string? ReportPath { get; init; }

    /// <summary>Default folder for "Export Thermo CSV..." (<c>--export-dir</c>), or null.</summary>
    public string? ExportDir { get; init; }

    /// <summary>Problems with the standalone options; empty when they parsed cleanly.</summary>
    public IReadOnlyList<string> Errors { get; init; } = Array.Empty<string>();

    public const string Usage = "SkylineCadenza.exe --report <report.parquet> [--export-dir <folder>]";

    public static LaunchOptions Parse(IReadOnlyList<string>? args)
    {
        args ??= Array.Empty<string>();
        if (args.Count == 0 || !args[0].StartsWith("--", StringComparison.Ordinal))
            return new LaunchOptions { SkylineArgs = args.ToArray() };

        string? report = null, exportDir = null;
        var errors = new List<string>();
        for (int i = 0; i < args.Count; i++)
        {
            string arg = args[i];
            // Accept both "--report path" and "--report=path".
            string name = arg;
            string? value = null;
            int eq = arg.IndexOf('=');
            if (arg.StartsWith("--", StringComparison.Ordinal) && eq > 0)
            {
                name = arg[..eq];
                value = arg[(eq + 1)..];
            }
            switch (name.ToLowerInvariant())
            {
                case "--report":
                case "--export-dir":
                    if (value is null)
                    {
                        if (i + 1 >= args.Count || args[i + 1].StartsWith("--", StringComparison.Ordinal))
                        {
                            errors.Add($"{name} needs a value.");
                            continue;
                        }
                        value = args[++i];
                    }
                    if (string.IsNullOrWhiteSpace(value))
                        errors.Add($"{name} needs a value.");
                    else if (name.Equals("--report", StringComparison.OrdinalIgnoreCase))
                        report = value;
                    else
                        exportDir = value;
                    break;
                default:
                    errors.Add($"Unknown option '{arg}'.");
                    break;
            }
        }
        return new LaunchOptions
        {
            Standalone = true,
            ReportPath = report,
            ExportDir = exportDir,
            Errors = errors,
        };
    }
}
