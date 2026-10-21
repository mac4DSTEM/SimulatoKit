"""SimulatoKit slice 4: any grid size, Potential(sampling=), beam tilt (abTEM 1.0.10, float64, numpy).

Writes Tests/SimulatoMultisliceTests/Fixtures/slice-4-pins.json and slice-4-arrays.bin: the grid abTEM derives from a
target sampling, and the exit waves of a plane wave and a probe through Al on grids no vDSP length matches (90 x 75),
untilted and tilted. abTEM states no tolerance for these paths; the tests use assert_allclose's defaults for float64
(rtol 1e-7) with atol 1e-12 of the wave's maximum.

Run:  refs/abtem/bin/python -I -B SimulatoKit/tools/slice-4/pins.py
"""

from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import abtem
from ase.build import bulk

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "Tests", "SimulatoMultisliceTests", "Fixtures")
BLOB = bytearray()
INDEX: dict = {}


def put(name, a):
    a = np.ascontiguousarray(a)
    raw = np.stack([a.real, a.imag], axis=-1).astype("<f8").tobytes()
    INDEX[name] = {"offset": len(BLOB), "count": int(a.size), "shape": list(a.shape), "kind": "c16"}
    BLOB.extend(raw)


abtem.config.set({"device": "cpu", "precision": "float64"})
atoms = bulk("Al", "fcc", a=4.05, cubic=True) * (2, 2, 3)
box = [float(x) for x in atoms.cell.lengths()]
energy = 200e3
cases = {"box": box, "energy": energy, "sliceThickness": 2.025, "atoms": [
    [int(z), [float(v) for v in p]] for z, p in zip(atoms.numbers, atoms.positions)]}

# 1. the grid from a target sampling
grids = []
for sampling in [0.1, 0.05, 0.123, (0.09, 0.11), 0.2025, 0.3]:
    p = abtem.Potential(atoms, sampling=sampling, slice_thickness=2.025, projection="infinite")
    grids.append({"sampling": sampling if isinstance(sampling, float) else list(sampling), "gpts": list(p.gpts)})
cases["grids"] = grids

# 2. exit waves on a 90 x 75 grid, untilted and tilted
gpts = (90, 75)
pot = abtem.Potential(atoms, gpts=gpts, slice_thickness=2.025, projection="infinite")
runs = {}
for label, tilt in [("plane0", (0.0, 0.0)), ("planeT", (1.5, -0.8)), ("probe0", (0.0, 0.0)), ("probeT", (2.0, 1.0))]:
    if label.startswith("plane"):
        w = abtem.PlaneWave(energy=energy, normalize=False, tilt=tilt, gpts=gpts, extent=pot.extent)
    else:
        w = abtem.Probe(energy=energy, semiangle_cutoff=20.0, defocus=0.0, tilt=tilt, gpts=gpts, extent=pot.extent)
        w = w  # centred at the scan position below
    if label.startswith("probe"):
        scan = abtem.CustomScan(np.array([[2.0, 3.0]]))
        out = w.multislice(pot, scan=scan).compute()
    else:
        out = w.multislice(pot).compute()
    arr = np.asarray(out.array).reshape(gpts)
    put(label, arr)
    runs[label] = {"tilt": list(tilt)}
cases["runs"] = runs
cases["probePosition"] = [2.0, 3.0]
cases["gpts"] = list(gpts)
cases["arrays"] = {"file": "slice-4-arrays.bin", "sha256": hashlib.sha256(bytes(BLOB)).hexdigest(), "index": INDEX}
with open(os.path.join(OUT, "slice-4-arrays.bin"), "wb") as f:
    f.write(BLOB)
with open(os.path.join(OUT, "slice-4-pins.json"), "w") as f:
    json.dump(cases, f)
print("abtem", abtem.__version__, "grids", grids)
