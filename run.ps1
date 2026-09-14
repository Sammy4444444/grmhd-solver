# Powershell build & run script for Windows environment

$ErrorActionPreference = "Stop"

Write-Host "=== [1/4] Membuat direktori build ===" -ForegroundColor Cyan
if (!(Test-Path -Path "build")) {
    New-Item -ItemType Directory -Path "build" | Out-Null
}

Set-Location -Path "build"

Write-Host "=== [2/4] Konfigurasi dan kompilasi CMake ===" -ForegroundColor Cyan
cmake ..
cmake --build . --config Release

Write-Host "=== [3/4] Menjalankan GRMHD Solver (C++) ===" -ForegroundColor Cyan
if (Test-Path "Release\grmhd_solver.exe") {
    .\Release\grmhd_solver.exe
} elseif (Test-Path "grmhd_solver.exe") {
    .\grmhd_solver.exe
} elseif (Test-Path "Debug\grmhd_solver.exe") {
    .\Debug\grmhd_solver.exe
} else {
    Write-Error "Executable grmhd_solver.exe tidak ditemukan."
}

Set-Location -Path ".."

Write-Host "=== [4/4] Menjalankan visualisasi Python & pengujian Julia ===" -ForegroundColor Cyan
if (Get-Command "python" -ErrorAction SilentlyContinue) {
    python scripts/plot_results.py output_2d.csv
} elseif (Get-Command "python3" -ErrorAction SilentlyContinue) {
    python3 scripts/plot_results.py output_2d.csv
} else {
    Write-Host "[INFO] Python tidak ditemukan di PATH." -ForegroundColor Yellow
}

if (Get-Command "julia" -ErrorAction SilentlyContinue) {
    Write-Host "--- Menjalankan prototype Julia ---" -ForegroundColor Green
    julia julia/test_metric.jl
} else {
    Write-Host "[INFO] Julia tidak ditemukan di PATH." -ForegroundColor Yellow
}

Write-Host "=== Selesai! Kerangka proyek GRMHD solver siap digunakan. ===" -ForegroundColor Green
