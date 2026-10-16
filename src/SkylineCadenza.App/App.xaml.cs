using System.Windows;
using SkylineCadenza.Core.SkylineRpc;

namespace SkylineCadenza.App;

public partial class App : Application
{
    /// <summary>
    /// Command-line args captured at startup. When Skyline launches the
    /// tool, <c>args[0]</c> is the named-pipe name passed via the
    /// <c>Arguments=$(SkylineConnection)</c> entry in
    /// <c>tool-inf/SkylineCadenza.properties</c>. The view-model uses it
    /// to attach to the running document.
    /// </summary>
    public static string[] LaunchArgs { get; private set; } = Array.Empty<string>();

    /// <summary>
    /// <see cref="LaunchArgs"/> parsed. When the first argument starts with
    /// <c>--</c> the app runs standalone (<c>--report</c> / <c>--export-dir</c>,
    /// no Skyline connection); otherwise it is the Skyline launch above.
    /// </summary>
    public static LaunchOptions Launch { get; private set; } = LaunchOptions.Parse(null);

    protected override void OnStartup(StartupEventArgs e)
    {
        LaunchArgs = e.Args ?? Array.Empty<string>();
        Launch = LaunchOptions.Parse(LaunchArgs);
        base.OnStartup(e);
    }
}
