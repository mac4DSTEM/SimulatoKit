#!/bin/zsh
# The benchmark table of docs/benchmarks.md: the kit side (swift test -c release) under tools/memcap.sh, then the abTEM
# side; one benchmark at a time machine-wide. Usage: tools/benchmarks/run.sh 128,256   (1024 and above: the nightly review)
# Output: tools/benchmarks/out/{kit,abtem}.txt, then `python3 tools/benchmarks/merge.py` writes the table.
set -eu
cd "$(dirname "$0")/../.."
sizes="${1:-128,256}"
mkdir -p tools/benchmarks/out
SIMULATO_BENCH_TABLE="$sizes" ../tools/memcap.sh -- swift test -c release --filter BenchmarkTableTests | grep '^BENCH|' > tools/benchmarks/out/kit.txt
../tools/memcap.sh -- ../refs/abtem/bin/python -I -B tools/benchmarks/abtem_bench.py "$sizes" | grep '^BENCH|' > tools/benchmarks/out/abtem.txt
