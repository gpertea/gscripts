param(
    [Parameter(Mandatory = $true)]
    [string]$Title
)

$ErrorActionPreference = 'Stop'

Add-Type -TypeDefinition @'
using System;
using System.Text;
using System.Runtime.InteropServices;

public static class VncWindowLookup {
    public delegate bool EnumWindowsProc(IntPtr hwnd, IntPtr parameter);

    [DllImport("user32.dll")]
    public static extern bool EnumWindows(EnumWindowsProc callback, IntPtr parameter);

    [DllImport("user32.dll")]
    public static extern uint GetWindowThreadProcessId(IntPtr hwnd, out uint processId);

    [DllImport("user32.dll", CharSet = CharSet.Unicode)]
    public static extern int GetWindowText(IntPtr hwnd, StringBuilder text, int maximum);

    [DllImport("user32.dll")]
    public static extern bool IsIconic(IntPtr hwnd);

    [DllImport("user32.dll")]
    public static extern bool IsWindowVisible(IntPtr hwnd);

    [DllImport("user32.dll")]
    public static extern bool ShowWindowAsync(IntPtr hwnd, int command);

    [DllImport("user32.dll")]
    public static extern bool SetForegroundWindow(IntPtr hwnd);
}
'@

$matches = [System.Collections.Generic.List[object]]::new()
$callback = [VncWindowLookup+EnumWindowsProc]{
    param($windowHandle, $parameter)

    $text = [Text.StringBuilder]::new(512)
    [void][VncWindowLookup]::GetWindowText($windowHandle, $text, $text.Capacity)
    if ($text.ToString() -ne $Title) {
        return $true
    }

    [uint32]$windowProcessId = 0
    [void][VncWindowLookup]::GetWindowThreadProcessId($windowHandle, [ref]$windowProcessId)
    $process = Get-Process -Id $windowProcessId -ErrorAction SilentlyContinue
    if ($null -eq $process -or $process.ProcessName -ne 'tigervncviewer') {
        return $true
    }

    $matches.Add([pscustomobject]@{
        Handle = $windowHandle
        Process = $process
        Minimized = [VncWindowLookup]::IsIconic($windowHandle)
        Visible = [VncWindowLookup]::IsWindowVisible($windowHandle)
    })
    return $true
}

[VncWindowLookup]::EnumWindows($callback, [IntPtr]::Zero) | Out-Null
if ($matches.Count -eq 0) {
    exit 1
}

$selected = $matches |
    Sort-Object @{ Expression = 'Visible'; Descending = $true },
                @{ Expression = 'Minimized'; Descending = $false },
                @{ Expression = { $_.Process.StartTime }; Descending = $true } |
    Select-Object -First 1

[void][VncWindowLookup]::ShowWindowAsync($selected.Handle, 9)
[void][VncWindowLookup]::SetForegroundWindow($selected.Handle)

Write-Output "Reusing existing TigerVNC viewer PID $($selected.Process.Id): $Title"
if ($matches.Count -gt 1) {
    Write-Warning "$($matches.Count) viewer windows match this session; run vnc-fix to diagnose stale clients."
}
exit 0
