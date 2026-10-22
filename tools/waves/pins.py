"""SimulatoKit waves detector: WavesDetector(gpts) over a scan (abTEM 1.0.10, float64, numpy).

Writes Tests/SimulatoKitParityTests/Fixtures/waves-pins.json and waves-arrays.bin: the exit waves of a 3 x 2 scan of a
probe through a small Al cell, as detected by WavesDetector() and by WavesDetector(gpts=...) for a crop on both axes, a
pad on both axes and a crop on one axis with a pad on the other (odd and even targets). abTEM states no tolerance for
`fft_interpolate`; the tests use assert_allclose's defaults for float64 (rtol 1e-7) with atol 1e-12 of the wave maximum.

Run:  References/python/abtem/bin/python -I -B SimulatoKit/tools/waves/pins.py
"""

from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import abtem
from ase.build import bulk

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "Tests", "SimulatoKitParityTests", "Fixtures")
BLOB = bytearray()
INDEX: dict = {}


def put(name, a):
    a = np.ascontiguousarray(a)
    raw = np.stack([a.real, a.imag], axis=-1).astype("<f8").tobytes()
    INDEX[name] = {"offset": len(BLOB), "count": int(a.size), "shape": list(a.shape), "kind": "c16"}
    BLOB.extend(raw)


assert abtem.__version__ == "1.0.10", abtem.__version__
abtem.config.set({"device": "cpu", "precision": "float64"})
atoms = bulk("Al", "fcc", a=4.05, cubic=True) * (1, 1, 2)
box = [float(x) for x in atoms.cell.lengths()]
gpts = 48
energy, semiangle, thickness = 200e3, 20.0, 2.025
pot = abtem.Potential(atoms, gpts=gpts, slice_thickness=thickness, projection="infinite")
probe = abtem.Probe(energy=energy, semiangle_cutoff=semiangle, gpts=gpts, extent=pot.extent)
scan = abtem.GridScan((0.0, 0.0), (4.05, 4.05), gpts=(3, 2), endpoint=False)

targets = [None, (33, 40), (81, 96), (30, 71), (64, 25), (10, 11)]
runs = []
for target in targets:
    if target is None:
        out = probe.scan(pot, scan, abtem.WavesDetector()).compute()
    else:
        # abTEM 1.0.10's `probe.scan(..., WavesDetector(gpts=))` fails (the detector does not declare its output shape, so the
        # preallocated measurement keeps the simulation grid: "could not broadcast (3,2,33,40) into (3,2,48,48)"). The
        # public `WavesDetector(gpts).detect(exit_waves)` on the scan's exit waves is the same code path per wave.
        exits = probe.multislice(pot, scan=scan).compute()
        out = abtem.WavesDetector(gpts=target).detect(exits).compute()
    arr = np.asarray(out.array)
    name = "full" if target is None else "gpts_%d_%d" % target
    put(name, arr)
    runs.append({"name": name, "gpts": None if target is None else list(target), "shape": list(arr.shape)})

cases = {
    "box": box, "gpts": gpts, "sliceThickness": thickness, "energy": energy, "semiangle": semiangle,
    "atoms": [[int(z), [float(v) for v in p]] for z, p in zip(atoms.numbers, atoms.positions)],
    "scan": {"start": [0.0, 0.0], "end": [4.05, 4.05], "gpts": [3, 2]},
    "runs": runs,
}
cases["arrays"] = {"file": "waves-arrays.bin", "sha256": hashlib.sha256(bytes(BLOB)).hexdigest(), "index": INDEX}
with open(os.path.join(OUT, "waves-arrays.bin"), "wb") as f:
    f.write(BLOB)
with open(os.path.join(OUT, "waves-pins.json"), "w") as f:
    json.dump(cases, f)
print("abtem", abtem.__version__, [(r["name"], r["shape"]) for r in runs])
