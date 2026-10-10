#!/bin/sh
# build.sh - build nightsky on Linux / macOS (POSIX).
#
# The Windows-only device-location code is guarded by #ifdef _WIN32, so on
# these platforms it compiles out and no Windows libraries are needed. Date/time
# handling uses standard C (time/gmtime), so it is identical across platforms.
# Location auto-detect is Windows-only; elsewhere nightsky falls back to
# nightsky.cfg (a "lat lon" line) then the built-in default.
#
# Requires: a C++20 compiler (g++ or clang++) and Python 3 (build-time only).
set -e

echo "[1/2] Generating binary catalogue nightsky.dat..."
python3 gencat.py

echo "[2/2] Compiling nightsky..."
CXX="${CXX:-c++}"
"$CXX" -std=c++20 -O2 -s -o nightsky nightsky.cpp

echo
echo "Build complete: ./nightsky (+ nightsky.dat)"
echo "Try:  ./nightsky 10"
