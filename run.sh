#!/usr/bin/env bash
# ==============================================================================
# Script otomatisasi kompilasi, eksekusi, dan visualisasi GRMHD Solver
# ==============================================================================

set -e

echo "=== [1/4] Membuat direktori build ==="
mkdir -p build
cd build

echo "=== [2/4] Konfigurasi dan kompilasi CMake ==="
cmake ..
cmake --build . --config Release

echo "=== [3/4] Menjalankan GRMHD Solver (C++) ==="
if [ -f "./grmhd_solver" ]; then
    ./grmhd_solver
elif [ -f "./Release/grmhd_solver.exe" ]; then
    ./Release/grmhd_solver.exe
elif [ -f "./grmhd_solver.exe" ]; then
    ./grmhd_solver.exe
elif [ -f "./Debug/grmhd_solver.exe" ]; then
    ./Debug/grmhd_solver.exe
else
    echo "[ERROR] Executable grmhd_solver tidak ditemukan."
    exit 1
fi

cd ..

echo "=== [4/4] Menjalankan visualisasi Python & pengujian Julia ==="
if command -v python3 &> /dev/null; then
    python3 scripts/plot_results.py output_2d.csv
elif command -v python &> /dev/null; then
    python scripts/plot_results.py output_2d.csv
else
    echo "[INFO] Python tidak terdeteksi. Silakan jalankan `python scripts/plot_results.py` secara manual."
fi

if command -v julia &> /dev/null; then
    echo "--- Menjalankan prototype Julia ---"
    julia julia/test_metric.jl
else
    echo "[INFO] Julia tidak terdeteksi. Anda dapat menguji file julia/test_metric.jl jika Julia sudah terinstal."
fi

echo "=== Selesai! Kerangka proyek GRMHD solver siap digunakan. ==="
