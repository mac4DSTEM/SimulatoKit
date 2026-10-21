"""SimulatoKit slice 4, scan side: GridScan(sampling=), centre of mass, Gaussian source size (abTEM 1.0.10, float64).

Writes Tests/SimulatoDetectorTests/Fixtures/slice-4-scan-pins.json and slice-4-scan-arrays.bin. abTEM's centre of mass
is complex64 (its own float32 accumulation), so the test bound there is 1e-5 of the largest value; the Gaussian filter
is scipy's in float64 (assert_allclose default, rtol 1e-7).

Run:  refs/abtem/bin/python -I -B SimulatoKit/tools/slice-4/scan_pins.py
"""

from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import abtem
from ase.build import bulk

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "Tests", "SimulatoDetectorTests", "Fixtures")
BLOB = bytearray()
INDEX: dict = {}


def put(name, a):
    a = np.ascontiguousarray(a)
    if np.iscomplexobj(a):
        raw, kind = np.stack([a.real, a.imag], axis=-1).astype("<f8").tobytes(), "c16"
    else:
        raw, kind = a.astype("<f8").tobytes(), "f8"
    INDEX[name] = {"offset": len(BLOB), "count": int(a.size), "shape": list(a.shape), "kind": kind}
    BLOB.extend(raw)


abtem.config.set({"device": "cpu", "precision": "float64"})
atoms = bulk("Al", "fcc", a=4.05, cubic=True) * (1, 1, 2)
box = [float(x) for x in atoms.cell.lengths()]
gpts = 64
pot = abtem.Potential(atoms, gpts=gpts, slice_thickness=2.025, projection="infinite")
probe = abtem.Probe(energy=200e3, semiangle_cutoff=20.0, gpts=gpts, extent=pot.extent)

# GridScan(sampling=): gpts from a sampling
scans = []
for start, end, sampling in [((0, 0), (4.05, 4.05), 0.5), ((0, 0), (4.05, 2.0), (0.3, 0.7)), ((1, 1), (3, 3), 0.25)]:
    s = abtem.GridScan(start, end, sampling=sampling)
    scans.append({"start": list(start), "end": list(end),
                  "sampling": sampling if not isinstance(sampling, tuple) else list(sampling), "gpts": list(s.gpts)})

# a 6 x 5 scan, pixelated patterns, centre of mass, Gaussian source size
scan = abtem.GridScan((0.0, 0.0), (4.05, 4.05), gpts=(6, 5), endpoint=False)
dp = probe.scan(pot, scan, abtem.PixelatedDetector(max_angle="cutoff")).compute()
put("dp", np.asarray(dp.array))
put("com_mrad", np.asarray(dp.center_of_mass(units="mrad").array))
put("com_inv", np.asarray(dp.center_of_mass(units="1/Å").array))
sigma = (0.8, 0.5)
put("dp_gauss", np.asarray(dp.gaussian_source_size(sigma).array))
sigma_big = 1.6
put("dp_gauss_big", np.asarray(dp.gaussian_source_size(sigma_big).array))
cases = {
    "box": box, "gpts": gpts, "sliceThickness": 2.025, "energy": 200e3, "semiangle": 20.0,
    "atoms": [[int(z), [float(v) for v in p]] for z, p in zip(atoms.numbers, atoms.positions)],
    "scans": scans, "scan": {"start": [0.0, 0.0], "end": [4.05, 4.05], "gpts": [6, 5]},
    "sigma": list(sigma), "sigmaBig": sigma_big, "scanSampling": list(scan.sampling),
    "dpShape": list(dp.array.shape),
}
cases["arrays"] = {"file": "slice-4-scan-arrays.bin", "sha256": hashlib.sha256(bytes(BLOB)).hexdigest(), "index": INDEX}
with open(os.path.join(OUT, "slice-4-scan-arrays.bin"), "wb") as f:
    f.write(BLOB)
with open(os.path.join(OUT, "slice-4-scan-pins.json"), "w") as f:
    json.dump(cases, f)
print("ok", scans, dp.array.shape, scan.sampling)
