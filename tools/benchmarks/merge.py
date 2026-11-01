"""Writes docs/benchmarks.md from tools/benchmarks/out/kit.txt and abtem.txt (run tools/benchmarks/run.sh first).

Run with the Python that has abTEM (its NumPy and SciPy versions go into the header):  .venv/bin/python tools/benchmarks/merge.py
"""

from __future__ import annotations

import os
import platform
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def read(name):
    rows = {}
    for line in open(os.path.join(HERE, "out", name)):
        _, function, size, backend, median, lo, hi = line.strip().split("|")
        rows[(function, int(size), backend)] = (float(median), float(lo), float(hi))
    return rows


def ms(t):
    m = t[0] * 1e3
    return f"{m:.3g} ({t[1] * 1e3:.3g}-{t[2] * 1e3:.3g})"


def sh(*cmd):
    return subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT).stdout.strip().splitlines()[0]


kit, ab = read("kit.txt"), read("abtem.txt")
functions = []
for f, _, _ in list(kit) + list(ab):
    if f not in functions:
        functions.append(f)
sizes = sorted({s for _, s, _ in kit})
out = []
out.append("# Benchmarks against abTEM\n")
out.append("Per-function timings of SimulatoKit against abTEM 1.0.10 on the same inputs (Al [001], 200 keV, a plane wave or a")
out.append("20 mrad probe, 10 slices of 2.025 Å, about 0.1 Å per pixel, infinite projection). One warm-up, then five timed runs;")
out.append("each cell is the median in milliseconds with (min-max). Written by `tools/benchmarks/run.sh` and `merge.py`.\n")
out.append(f"- Machine: {sh('sysctl', '-n', 'machdep.cpu.brand_string')}, macOS {sh('sw_vers', '-productVersion')}, "
           f"{sh('swift', '--version')}")
py = subprocess.run([os.environ.get("ABTEM_PYTHON", sys.executable), "-I", "-c",
                    "import sys, numpy, scipy; print(sys.version.split()[0], numpy.__version__, scipy.__version__)"],
                   capture_output=True, text=True).stdout.split()
out.append(f"- Python {py[0]} (NumPy {py[1]}, SciPy {py[2]}), abTEM 1.0.10 (CPU, NumPy's FFT; no GPU)")
out.append("- Kit CPU is float64 on Accelerate (one thread per call); Metal is float32 on MPSGraph; abTEM is shown in float64 and float32")
load = {}
for line in open(os.path.join(HERE, "out", "load.txt")):  # run.sh: "start { 1m 5m 15m }" and "end { ... }"
    key, _, rest = line.partition(" ")
    load[key] = rest.strip().strip("{}").strip()
out.append(f"- Load average (1, 5, 15 min): {load.get('start', '?')} at the start, {load.get('end', '?')} at the end; "
           "the spread columns show the noise; treat ratios within 20% as equal")
out.append(f"- Sizes measured: {', '.join(f'{s}²' for s in sizes)}\n")
out.append("Re-run from a fresh clone of the kit, on a quiet Mac (`MEMORY_CAP` may name a memory-cap command prefix instead of `none`):\n")
out.append("```sh")
out.append("python3.12 -m venv .venv && .venv/bin/pip install -r tools/requirements.txt")
out.append("MEMORY_CAP=none tools/benchmarks/run.sh 64,128,192,256,512")
out.append(".venv/bin/python tools/benchmarks/merge.py   # writes this file")
out.append("```\n")

for f in functions:
    out.append(f"## {f}\n")
    out.append("| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |")
    out.append("|---|---|---|---|---|---|---|")
    fsizes = sorted({s for (g, s, _) in list(kit) + list(ab) if g == f})
    for s in fsizes:
        c, m = kit.get((f, s, "cpu")), kit.get((f, s, "metal"))
        a64, a32 = ab.get((f, s, "float64")), ab.get((f, s, "float32"))
        best = min([x for x in (c, m) if x], key=lambda t: t[0], default=None)
        r64 = f"{a64[0] / c[0]:.2f}x" if a64 and c else ""
        r32 = f"{a32[0] / best[0]:.2f}x" if a32 and best else ""
        out.append(f"| {s}² | {ms(c) if c else ''} | {ms(m) if m else ''} | {ms(a64) if a64 else ''} | {ms(a32) if a32 else ''} | {r64} | {r32} |")
    out.append("\n(ratios are abTEM time / kit time: above 1 the kit is faster)\n")

# backend decision
out.append("## Backend decision per function\n")
out.append("The smallest measured size at which Metal beats the CPU path. The break-even is bracketed by the sizes measured, not interpolated.\n")
out.append("| function | decision | measured |")
out.append("|---|---|---|")
by = defaultdict(dict)
for (f, s, b), t in kit.items():
    by[f].setdefault(s, {})[b] = t[0]
for f, per in by.items():
    both = {s: v for s, v in per.items() if "cpu" in v and "metal" in v}
    if not both and any("metal" in v for v in per.values()):
        why = "Metal only (the GPU build of the Metal row above; its CPU counterpart is potential build + multislice)"
        out.append(f"| {f} | {why} | |")
        continue
    if not both:
        why = "CPU only: no Metal path yet (see the issues if it is a work item), or the work is CPU-side by nature"
        out.append(f"| {f} | {why} | |")
        continue
    wins = [s for s in sorted(both) if both[s]["metal"] < both[s]["cpu"] * 0.9]
    detail = ", ".join(f"{s}²: cpu {both[s]['cpu'] * 1e3:.3g} ms / metal {both[s]['metal'] * 1e3:.3g} ms" for s in sorted(both))
    if wins and all(both[s]["metal"] < both[s]["cpu"] * 0.9 for s in sorted(both) if s >= wins[0]):
        prev = [s for s in sorted(both) if s < wins[0]]
        out.append(f"| {f} | Metal above {wins[0]}² (break-even between {prev[-1] if prev else '<' + str(wins[0])}² and {wins[0]}²) | {detail} |")
    else:
        out.append(f"| {f} | CPU only at the sizes measured | {detail} |")
with open(os.path.join(ROOT, "docs", "benchmarks.md"), "w") as f:
    f.write("\n".join(out) + "\n")
print("docs/benchmarks.md written")
