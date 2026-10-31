"""SimulatoKit line-scan defaults: LineScan(end=None) and LineScan.match_probe (abTEM 1.0.10, numpy).

Prints the numbers that Tests/SimulatoProbeTests/LineScanDefaultsTests.swift holds inline (no fixture file): line scans
without an end point or a sampling, matched to a probe (a 20 mrad probe at 200 keV and a 30 mrad probe at 80 keV on a
12.3 x 9.7 A grid). match_probe puts a missing end at (0, extent[1]) of the probe and, without a sampling, sets 0.99 x the
Nyquist sampling of the aperture, which recomputes the number of positions. The test holds the end points and counts
exact and the sampling at rtol 1e-7 (assert_allclose defaults; abTEM states no tolerance).

Run from the repo root:  .venv/bin/python -I -B tools/pins/linescan-defaults/pins.py
"""

from abtem import LineScan, Probe

probes = {
    "p20": Probe(semiangle_cutoff=20.0, energy=200e3, extent=(12.3, 9.7), gpts=(64, 48)),
    "p30": Probe(semiangle_cutoff=30.0, energy=80e3, extent=(12.3, 9.7), gpts=(64, 48)),
}

cases = {
    "default": dict(),
    "start only": dict(start=(1.0, 2.0)),
    "no end, no endpoint": dict(endpoint=False),
    "end, no sampling": dict(start=(1.0, 1.0), end=(8.0, 5.0)),
    "end and gpts": dict(end=(8.0, 5.0), gpts=17),
    "no end, sampling": dict(sampling=0.3),
    "no end, gpts": dict(gpts=10),
}

for pname, probe in probes.items():
    for cname, kwargs in cases.items():
        scan = LineScan(**kwargs)
        scan.match_probe(probe)
        print(pname, "|", cname, "| start", scan.start, "end", scan.end, "gpts", scan.gpts, "sampling", repr(scan.sampling))
