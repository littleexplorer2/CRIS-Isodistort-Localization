$ErrorActionPreference = "Stop"

$crisRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$venvRoot = Join-Path $crisRoot ".venv"
$venvConfig = Join-Path $venvRoot "pyvenv.cfg"
$sitePackages = Join-Path $venvRoot "Lib\site-packages"
$venvScripts = Join-Path $venvRoot "Scripts"

if (-not (Test-Path -LiteralPath $venvConfig -PathType Leaf)) {
    throw "CRIS virtual environment is missing: $venvConfig"
}
if (-not (Test-Path -LiteralPath $sitePackages -PathType Container)) {
    throw "CRIS site-packages directory is missing: $sitePackages"
}

$basePythonLine = Get-Content -LiteralPath $venvConfig |
    Where-Object { $_ -match '^executable\s*=\s*(.+)$' } |
    Select-Object -First 1
if (-not $basePythonLine) {
    throw "No base Python executable is recorded in $venvConfig"
}
$basePython = ($basePythonLine -split '=', 2)[1].Trim()
if (-not (Test-Path -LiteralPath $basePython -PathType Leaf)) {
    throw "The base Python recorded by the CRIS environment is missing: $basePython"
}

# On this OneDrive installation, a Python executable whose image is inside the
# CRIS directory receives Wsl/E_ACCESSDENIED. Start the recorded base Python
# outside OneDrive, while loading every third-party package from CRIS/.venv.
$env:VIRTUAL_ENV = $venvRoot
$env:PYTHONNOUSERSITE = "1"
$pythonPathEntries = @($sitePackages, (Join-Path $crisRoot "ISODISTORT"))
$env:PYTHONPATH = $pythonPathEntries -join [IO.Path]::PathSeparator
$env:PATH = $venvScripts + [IO.Path]::PathSeparator + $env:PATH

if ($args.Count -eq 0) {
    $pythonArguments = @((Join-Path $crisRoot "ISODISTORT\scripts\main_web.py"))
} else {
    $pythonArguments = @($args)
}

# -S prevents the external base installation from adding its global
# site-packages. PYTHONPATH above supplies the physical CRIS environment.
& $basePython -S @pythonArguments
exit $LASTEXITCODE
