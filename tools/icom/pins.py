"""SimulatoKit iCOM slice: Images.integrate_gradient, integrated centre of mass, integrate_disc, Images.tile
(abTEM 1.0.10, float64, numpy).

Writes Tests/SimulatoKitParityTests/Fixtures/icom-pins.json and icom-arrays.bin:
  1. planted gradient fields: the analytic gradient of known smooth periodic potentials on even, odd and non-square
     grids with unequal sampling, and abTEM's `Images.integrate_gradient` of each (it recovers the potential up to a
     constant, then abTEM subtracts the minimum); plus a random gradient (all frequencies, Nyquist included);
  2. a small simulated scan (12 x 10, Al, pixelated patterns): abTEM's centre of mass in 1/Å and mrad, and its
     `integrated_center_of_mass` of them; the kit integrates abTEM's centre of mass (strict) and its own (looser);
  3. `integrate_disc` of a random image at positions that wrap, with sum and mean, and `Images.tile`.
abTEM states no tolerance for these functions; the tests use assert_allclose's defaults for float64 (rtol 1e-7) with
atol 1e-12 of the maximum where a value is zero in one arithmetic.

Run:  References/python/abtem/bin/python -I -B SimulatoKit/tools/icom/pins.py
"""

from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import abtem
from abtem.measurements import Images, integrate_disc
from ase.build import bulk

assert abtem.__version__ == "1.0.10", abtem.__version__

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "Tests", "SimulatoKitParityTests", "Fixtures")
BLOB = bytearray()
INDEX: dict = {}


def put(name, a):
    a = np.ascontiguousarray(a)
    if np.iscomplexobj(a):
        a = a.astype(np.complex128)
        raw, kind = np.stack([a.real, a.imag], axis=-1).astype("<f8").tobytes(), "c16"
    else:
        raw, kind = a.astype("<f8").tobytes(), "f8"
    INDEX[name] = {"offset": len(BLOB), "count": int(a.size), "shape": list(a.shape), "kind": kind}
    BLOB.extend(raw)


abtem.config.set({"device": "cpu", "precision": "float64"})
rng = np.random.default_rng(20261022)
cases: dict = {}

# 1. planted gradients
grads = []
for label, (nx, ny), (sx, sy) in [("even", (16, 12), (0.3, 0.5)), ("odd", (15, 11), (0.4, 0.4)),
                                  ("mixed", (20, 9), (0.25, 0.6)), ("tiny", (2, 3), (1.0, 1.0))]:
    ex, ey = nx * sx, ny * sy
    x = np.arange(nx) * sx
    y = np.arange(ny) * sy
    X, Y = np.meshgrid(x, y, indexing="ij")
    T = np.zeros((nx, ny))
    gx = np.zeros((nx, ny))
    gy = np.zeros((nx, ny))
    # components (m, n) inside the grid's non-Nyquist band, amplitude and phase from the generator
    for m, n in [(1, 0), (0, 1), (1, 1), (2, -1), (1, 2)]:
        if abs(m) * 2 >= nx or abs(n) * 2 >= ny:
            continue
        a, ph = rng.uniform(0.5, 1.5), rng.uniform(0, 2 * np.pi)
        kx, ky = m / ex, n / ey
        arg = 2 * np.pi * (kx * X + ky * Y) + ph
        T += a * np.cos(arg)
        gx += -a * 2 * np.pi * kx * np.sin(arg)
        gy += -a * 2 * np.pi * ky * np.sin(arg)
    g = gx + 1j * gy
    out = np.asarray(Images(g, sampling=(sx, sy)).integrate_gradient().array)
    put(f"grad_{label}_T", T)
    put(f"grad_{label}_g", g)
    put(f"grad_{label}_out", out)
    grads.append({"label": label, "shape": [nx, ny], "sampling": [sx, sy]})
# random gradient: all frequencies including Nyquist
for label, (nx, ny), (sx, sy) in [("rand", (14, 10), (0.2, 0.35)), ("randOdd", (13, 9), (0.5, 0.2))]:
    g = rng.normal(size=(nx, ny)) + 1j * rng.normal(size=(nx, ny))
    out = np.asarray(Images(g, sampling=(sx, sy)).integrate_gradient().array)
    put(f"grad_{label}_g", g)
    put(f"grad_{label}_out", out)
    grads.append({"label": label, "shape": [nx, ny], "sampling": [sx, sy]})
cases["gradients"] = grads

# 2. a small simulated scan
atoms = bulk("Al", "fcc", a=4.05, cubic=True) * (1, 1, 2)
gpts = 64
pot = abtem.Potential(atoms, gpts=gpts, slice_thickness=2.025, projection="infinite")
probe = abtem.Probe(energy=200e3, semiangle_cutoff=20.0, gpts=gpts, extent=pot.extent)
scan = abtem.GridScan((0.0, 0.0), (4.05, 4.05), gpts=(12, 10), endpoint=False)
dp = probe.scan(pot, scan, abtem.PixelatedDetector(max_angle="cutoff")).compute()
scan_sampling = [float(s) for s in scan.sampling]
for units, tag in [("1/Å", "inv"), ("mrad", "mrad")]:
    com = dp.center_of_mass(units=units)
    put(f"scan_com_{tag}", np.asarray(com.array))
    put(f"scan_icom_{tag}", np.asarray(com.integrate_gradient().array))
    cases[f"comDtype_{tag}"] = str(np.asarray(com.array).dtype)
put("scan_icom_chain", np.asarray(dp.integrated_center_of_mass().array))
cases["scan"] = {
    "box": [float(x) for x in atoms.cell.lengths()], "gpts": gpts, "sliceThickness": 2.025, "energy": 200e3,
    "semiangle": 20.0, "start": [0.0, 0.0], "end": [4.05, 4.05], "scanGpts": [12, 10], "sampling": scan_sampling,
    "atoms": [[int(z), [float(v) for v in p]] for z, p in zip(atoms.numbers, atoms.positions)],
}

# 3. integrate_disc and tile
nx, ny, sx, sy = 17, 13, 0.21, 0.33
img = rng.normal(size=(nx, ny))
im = Images(img, sampling=(sx, sy))
put("disc_image", img)
discs = []
for pos, radius in [((1.7, 2.1), 0.5), ((0.05, 0.1), 0.6), ((3.5, 4.2), 0.7), ((-0.2, 4.0), 0.45), ((3.4, 0.0), 1.1),
                    ((1.0, 1.0), 0.08), ((1.78, 2.14), 2.0)]:
    s = float(integrate_disc(im, pos, radius))
    m = float(integrate_disc(im, pos, radius, return_mean=True))
    discs.append({"position": list(pos), "radius": radius, "sum": s, "mean": m})
cases["disc"] = {"shape": [nx, ny], "sampling": [sx, sy], "cases": discs}
tiled = np.asarray(im.tile((2, 3)).array)
put("tile_2_3", tiled)
cases["tile"] = {"repetitions": [2, 3], "shape": list(tiled.shape)}

cases["arrays"] = {"file": "icom-arrays.bin", "sha256": hashlib.sha256(bytes(BLOB)).hexdigest(), "index": INDEX}
with open(os.path.join(OUT, "icom-arrays.bin"), "wb") as f:
    f.write(BLOB)
with open(os.path.join(OUT, "icom-pins.json"), "w") as f:
    json.dump(cases, f)
print("abtem", abtem.__version__, "com dtypes", cases["comDtype_inv"], cases["comDtype_mrad"], "sampling", scan_sampling)
