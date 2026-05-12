# Test coverage report for bromesh (Windows / MSVC only).
#
# Runs the bromesh_test suite (via ctest) under OpenCppCoverage and emits an
# HTML report at build/coverage/index.html.
#
# IMPORTANT: bromesh tests MUST run in Release. meshoptimizer's Debug
# assertions trigger a modal abort() dialog on Windows that hangs the run.
#
# Requirements:
#   - OpenCppCoverage installed (winget install OpenCppCoverage.OpenCppCoverage)
#   - Release build present at build/tests/Release/bromesh_test.exe (PDBs needed)
#     Reconfigure with: cmake -B build -DCMAKE_MSVC_DEBUG_INFORMATION_FORMAT=ProgramDatabase
#     Build with:       cmake --build build --config Release --target bromesh_test
#
# Usage:
#   pwsh scripts/coverage.ps1
#   pwsh scripts/coverage.ps1 -Output build/cov   # custom output dir

[CmdletBinding()]
param(
    [string]$Output = 'build/coverage'
)

$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

$occ = "C:\Program Files\OpenCppCoverage\OpenCppCoverage.exe"
if (-not (Test-Path $occ)) {
    throw "OpenCppCoverage not found at $occ. Install: winget install OpenCppCoverage.OpenCppCoverage"
}

$testExe = Join-Path $root 'build\tests\Release\bromesh_test.exe'
if (-not (Test-Path $testExe)) {
    throw "$testExe not found. Build first: cmake --build build --config Release --target bromesh_test"
}

$outAbs = if ([System.IO.Path]::IsPathRooted($Output)) { $Output } else { Join-Path $root $Output }
if (Test-Path $outAbs) { Remove-Item -Recurse -Force $outAbs }

& $occ `
    --sources "$root\src" `
    --modules "$testExe" `
    --export_type "html:$outAbs" `
    --working_dir $root `
    --quiet `
    -- $testExe

Write-Host ""
Write-Host "Coverage report: $outAbs\index.html"
