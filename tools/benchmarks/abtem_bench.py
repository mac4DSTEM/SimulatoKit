"""The Python side of docs/benchmarks.md: abTEM 1.0.10 (float64 and float32, CPU, numpy/scipy FFT) on the same inputs as
Tests/SimulatoDetectorTests/BenchmarkTableTests.swift. One warm-up, then five timed runs; prints
`BENCH|function|size|backend|median|min|max` in seconds.

Run:  ../refs/abtem/bin/python -I -B tools/benchmarks/abtem_bench.py 128,256
"""

from __future__ import annotations

import sys
import time

import numpy as np
import abtem
from ase import Atoms


def report(function, size, backend, body):
    body()
    t = []
    for _ in range(5):
        a = time.perf_counter()
        body()
        t.append(time.perf_counter() - a)
    t.sort()
    print(f"BENCH|{function}|{size}|{backend}|{t[2]:.6f}|{t[0]:.6f}|{t[4]:.6f}", flush=True)


def structure(n):
    cell = Atoms("Al4", positions=[(0, 0, 0), (0, 2.025, 2.025), (2.025, 0, 2.025), (2.025, 2.025, 0)], cell=(4.05,) * 3, pbc=True)
    cells = max(1, int(round(n * 0.1 / 4.05)))
    return cell * (cells, cells, 5)


def run(n, precision):
    abtem.config.set({"device": "cpu", "precision": precision, "dask.lazy": False})
    atoms = structure(n)
    energy = 200e3
    pot = abtem.Potential(atoms, gpts=n, slice_thickness=2.025, projection="infinite")
    built = pot.build(lazy=False)
    wave = abtem.PlaneWave(energy=energy, normalize=False, gpts=n, extent=pot.extent)
    for m in (n, n - 6):
        x = (np.sin(np.arange(m * m, dtype=np.float64)).reshape(m, m) + 0j).astype(np.complex128 if precision == "float64" else np.complex64)
        def fft():
            nonlocal x
            x = np.fft.ifft2(np.fft.fft2(x))
        report("fft2+ifft2", m, precision, fft)
    report("potential build (infinite)", n, precision, lambda: pot.build(lazy=False))
    report("multislice 10 slices (plane wave, built potential)", n, precision,
           lambda: wave.multislice(built, lazy=False))
    probe = abtem.Probe(energy=energy, semiangle_cutoff=20.0, gpts=n, extent=pot.extent)
    report("probe build", n, precision, lambda: probe.build((1.0, 1.0), lazy=False))
    scan = abtem.GridScan((0, 0), (4.05, 4.05), gpts=(8, 8), endpoint=False)
    det = abtem.AnnularDetector(inner=40, outer=90)
    report("STEM scan 8x8, annular", n, precision, lambda: probe.scan(built, scan, det, lazy=False))
    fp = abtem.FrozenPhonons(atoms, num_configs=4, sigmas=0.08, seed=1)
    fpot = abtem.Potential(fp, gpts=n, slice_thickness=2.025, projection="infinite")
    report("frozen phonons 4 x multislice, built per configuration", n, precision,
           lambda: wave.multislice(fpot, lazy=False))


if __name__ == "__main__":
    for n in [int(s) for s in sys.argv[1].split(",")]:
        for precision in ("float64", "float32"):
            run(n, precision)
