using SkylineCadenza.Core.SkylineRpc;
using Xunit;

namespace SkylineCadenza.Tests;

/// <summary>
/// The app's command line: a Skyline pipe name keeps today's launch; a first
/// argument starting with "--" means a standalone launch.
/// </summary>
public class LaunchOptionsTests
{
    [Fact]
    public void No_arguments_is_a_skyline_launch()
    {
        var o = LaunchOptions.Parse(Array.Empty<string>());
        Assert.False(o.Standalone);
        Assert.Empty(o.SkylineArgs);
        Assert.Null(o.ReportPath);
    }

    [Fact]
    public void Pipe_name_is_passed_through_untouched()
    {
        var args = new[] { "SkylineToolService-1234", "extra" };
        var o = LaunchOptions.Parse(args);
        Assert.False(o.Standalone);
        Assert.Equal(args, o.SkylineArgs);
        Assert.Null(o.ReportPath);
        Assert.Null(o.ExportDir);
        Assert.Empty(o.Errors);
    }

    [Fact]
    public void Report_and_export_dir_make_a_standalone_launch()
    {
        var o = LaunchOptions.Parse(new[]
        {
            "--report", @"C:\data\run 1\report.parquet", "--export-dir", @"C:\out\assay",
        });
        Assert.True(o.Standalone);
        Assert.Empty(o.SkylineArgs);
        Assert.Equal(@"C:\data\run 1\report.parquet", o.ReportPath);
        Assert.Equal(@"C:\out\assay", o.ExportDir);
        Assert.Empty(o.Errors);
    }

    [Fact]
    public void Equals_form_and_option_order_are_accepted()
    {
        var o = LaunchOptions.Parse(new[] { "--export-dir=C:\\out", "--report=r.parquet" });
        Assert.True(o.Standalone);
        Assert.Equal("r.parquet", o.ReportPath);
        Assert.Equal("C:\\out", o.ExportDir);
        Assert.Empty(o.Errors);
    }

    [Fact]
    public void Export_dir_alone_is_standalone_without_a_report()
    {
        var o = LaunchOptions.Parse(new[] { "--export-dir", "out" });
        Assert.True(o.Standalone);
        Assert.Null(o.ReportPath);
        Assert.Equal("out", o.ExportDir);
    }

    [Fact]
    public void Missing_value_and_unknown_option_are_reported()
    {
        var o = LaunchOptions.Parse(new[] { "--report", "--bogus" });
        Assert.True(o.Standalone);
        Assert.Null(o.ReportPath);
        Assert.Equal(2, o.Errors.Count);
        Assert.Contains(o.Errors, e => e.Contains("--report"));
        Assert.Contains(o.Errors, e => e.Contains("--bogus"));
    }
}
