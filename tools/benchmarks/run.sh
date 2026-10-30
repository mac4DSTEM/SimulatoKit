#!/bin/zsh
# The benchmark table of docs/benchmarks.md: the kit side (swift test -c release) under a memory cap, then the abTEM
# side; one benchmark at a time machine-wide. Usage: tools/benchmarks/run.sh 128,256   (1024 and above: the nightly review)
# Output: tools/benchmarks/out/{kit,abtem,load}.txt (ignored by git), then `.venv/bin/python tools/benchmarks/merge.py` writes the table.
set -eu
# MEMORY_CAP: the command prefix that caps the run's memory (a wrapper that stops the process tree past its limit, with its
# argument separator), as an absolute path; MEMORY_CAP=none runs without a cap.
# ABTEM_PYTHON: the Python with abTEM 1.0.10 installed (default: .venv/bin/python, set up as in tools/README.md).
: ${MEMORY_CAP:?set MEMORY_CAP to a memory-cap command prefix, or none}
cap=(); [[ $MEMORY_CAP == none ]] || cap=(${=MEMORY_CAP})
cd "$(dirname "$0")/../.."
sizes="${1:-128,256}"
abtem_python="${ABTEM_PYTHON:-.venv/bin/python}"
[[ -x $abtem_python ]] || { echo "run.sh: no Python at $abtem_python; set up .venv (tools/README.md) or set ABTEM_PYTHON" >&2; exit 2; }
mkdir -p tools/benchmarks/out
echo "start $(sysctl -n vm.loadavg)" > tools/benchmarks/out/load.txt
SIMULATO_BENCH_TABLE="$sizes" $cap swift test -c release --filter BenchmarkTableTests | grep '^BENCH|' > tools/benchmarks/out/kit.txt
$cap "$abtem_python" -I -B tools/benchmarks/abtem_bench.py "$sizes" | grep '^BENCH|' > tools/benchmarks/out/abtem.txt
echo "end $(sysctl -n vm.loadavg)" >> tools/benchmarks/out/load.txt
