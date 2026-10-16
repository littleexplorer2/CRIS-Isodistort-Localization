param([string]$Latexmk = 'latexmk')
$ErrorActionPreference = 'Stop'
$libraryRoot = Split-Path $PSScriptRoot -Parent
$cacheRoot = Join-Path $libraryRoot 'audit\.texcache'
$savedEnvironment = @{}
foreach ($key in @('TEXMFVAR', 'TEXMFCONFIG', 'TEXMFCACHE', 'LC_ALL')) {
    $savedEnvironment[$key] = [Environment]::GetEnvironmentVariable($key, 'Process')
}
try {
    $env:TEXMFVAR = Join-Path $cacheRoot 'var'
    $env:TEXMFCONFIG = Join-Path $cacheRoot 'config'
    $env:TEXMFCACHE = Join-Path $cacheRoot 'cache'
    $env:LC_ALL = 'C'
    Push-Location -LiteralPath $PSScriptRoot
    try {
        & $Latexmk -quiet -xelatex -interaction=nonstopmode -halt-on-error -outdir=build -latexoption=-no-shell-escape unsupported_algorithm_proofs.tex isodistort_theoretical_foundations.tex
        if ($LASTEXITCODE -ne 0) { throw "LaTeX build failed with exit code $LASTEXITCODE" }
        foreach ($name in @('unsupported_algorithm_proofs', 'isodistort_theoretical_foundations')) {
            Copy-Item -LiteralPath (Join-Path $PSScriptRoot "build\$name.pdf") -Destination (Join-Path $PSScriptRoot "$name.pdf")
        }
    } finally { Pop-Location }
} finally {
    foreach ($key in $savedEnvironment.Keys) {
        [Environment]::SetEnvironmentVariable($key, $savedEnvironment[$key], 'Process')
    }
}
