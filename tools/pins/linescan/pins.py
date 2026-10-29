"""SimulatoKit line scans and line profiles (pins): abTEM 1.0.10, float64.

Pins `LineScan` (points, sampling, endpoint, `add_margin`, `at_position`, fractional coordinates), a `STEM.scan` over a
LineScan with an annular detector, the `RealSpaceLineProfiles` width and interpolate (on the scan's profile and on planted
arrays, orders 0 to 5, both endpoints), `Images.interpolate_line` and `interpolate_line_at_position` (width, margin,
fractional, orders), and the same on a `DiffractionPatterns` (reciprocal-space profiles).

abTEM states no tolerance for these paths. The tests use `assert_allclose`'s float64 defaults (rtol 1e-7, with atol 1e-12 of the
maximum where a value is zero in one arithmetic); abTEM stores `interpolate` and `width` results as float32, which is
itself inside rtol 1e-7 of the float64 value to 6e-8. The scan runs through the kit's potential with abTEM's float32 pi
(`LobatoParameterization.PiConstant.abTEMFloat32`), as the other pins do.

Writes Tests/SimulatoKitParityTests/Fixtures/linescan-pins.json and linescan-arrays.bin (floats as repr in the JSON, arrays
as little-endian float64 in the sidecar, with its sha256).

Run from the repo root:  .venv/bin/python -I -B tools/pins/linescan/pins.py
"""

from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import scipy
import abtem
from abtem.core.axes import UnknownAxis
from ase.build import bulk

assert abtem.__version__ == "1.0.10", abtem.__version__
assert np.__version__ == "2.5.4", np.__version__
assert scipy.__version__ == "1.18.1", scipy.__version__

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "..", "Tests", "SimulatoKitParityTests", "Fixtures")
BLOB = bytearray()
INDEX: dict = {}


def put(name, a):
    a = np.ascontiguousarray(a, dtype="<f8")
    INDEX[name] = {"offset": len(BLOB), "count": int(a.size), "shape": list(a.shape), "kind": "f8"}
    BLOB.extend(a.tobytes())


abtem.config.set({"device": "cpu", "precision": "float64"})
rng = np.random.default_rng(20261022)

atoms = bulk("Al", "fcc", a=4.05, cubic=True) * (1, 1, 2)
box = [float(x) for x in atoms.cell.lengths()]
gpts = 64
pot = abtem.Potential(atoms, gpts=gpts, slice_thickness=2.025, projection="infinite")
cases: dict = {
    "versions": {"abtem": abtem.__version__, "numpy": np.__version__, "scipy": scipy.__version__},
    "box": box, "gpts": gpts, "sliceThickness": 2.025, "energy": 200e3, "semiangle": 20.0,
    "atoms": [[int(z), [float(v) for v in p]] for z, p in zip(atoms.numbers, atoms.positions)],
}

# ---- 1. LineScan: construction, margin, at_position, fractional ------------------------------------------------------
line_cases = []


def record(label, scan, **args):
    positions = scan.get_positions()
    put(f"line.{label}.positions", positions)
    line_cases.append({"label": label, "args": args, "start": list(scan.start), "end": list(scan.end),
                       "gpts": int(scan.gpts), "sampling": float(scan.sampling), "endpoint": bool(scan.endpoint),
                       "extent": float(scan.extent), "angle": float(scan.angle)})


record("gpts7", abtem.LineScan((0, 0), (4.05, 2.0), gpts=7), start=[0, 0], end=[4.05, 2.0], gpts=7)
record("sampling_end", abtem.LineScan((0, 0), (4.05, 2.0), sampling=0.5), start=[0, 0], end=[4.05, 2.0], sampling=0.5)
record("sampling_noend", abtem.LineScan((0, 0), (4.05, 2.0), sampling=0.5, endpoint=False),
       start=[0, 0], end=[4.05, 2.0], sampling=0.5, endpoint=False)
record("both", abtem.LineScan((0, 0), (4.05, 2.0), gpts=5, sampling=0.5), start=[0, 0], end=[4.05, 2.0], gpts=5, sampling=0.5)
record("gpts1", abtem.LineScan((0.3, 0.2), (4.05, 2.0), gpts=1), start=[0.3, 0.2], end=[4.05, 2.0], gpts=1)
record("gpts1_noend", abtem.LineScan((0.3, 0.2), (4.05, 2.0), gpts=1, endpoint=False),
       start=[0.3, 0.2], end=[4.05, 2.0], gpts=1, endpoint=False)
record("margin_gpts", abtem.LineScan((1, 1), (4.05, 2.0), gpts=9).add_margin(0.3),
       start=[1, 1], end=[4.05, 2.0], gpts=9, margin=[0.3, 0.3])
record("margin_asym", abtem.LineScan((0.5, 0.2), (3.0, 3.5), sampling=0.25).add_margin((0.2, 0.7)),
       start=[0.5, 0.2], end=[3.0, 3.5], sampling=0.25, margin=[0.2, 0.7])
record("margin_noend", abtem.LineScan((0.5, 0.2), (3.0, 3.5), gpts=12, endpoint=False).add_margin(0.45),
       start=[0.5, 0.2], end=[3.0, 3.5], gpts=12, endpoint=False, margin=[0.45, 0.45])
record("at_pos_gpts", abtem.LineScan.at_position((2.0, 2.0), extent=3.0, angle=30.0, gpts=11),
       center=[2.0, 2.0], extent=3.0, angle=30.0, gpts=11, atPosition=True)
record("at_pos_sampling", abtem.LineScan.at_position((1.5, 2.5), extent=2.2, angle=90.0, sampling=0.4, endpoint=False),
       center=[1.5, 2.5], extent=2.2, angle=90.0, sampling=0.4, endpoint=False, atPosition=True)
record("at_pos_default", abtem.LineScan.at_position((2.0, 2.0), gpts=4), center=[2.0, 2.0], gpts=4, atPosition=True)
record("fractional", abtem.LineScan((0.0, 0.5), (1.0, 0.25), gpts=10, fractional=True, potential=pot),
       start=[0.0, 0.5], end=[1.0, 0.25], gpts=10, fractional=True, extent2=[float(e) for e in pot.extent])
cases["lineScans"] = line_cases

# ---- 2. STEM.scan over a LineScan, an annular detector ---------------------------------------------------------------
probe = abtem.Probe(energy=200e3, semiangle_cutoff=20.0, gpts=gpts, extent=pot.extent)
scan_specs = []
for label, scan, spec in [
    ("scanA", abtem.LineScan((0.0, 0.0), (4.05, 2.0), gpts=9), {"start": [0, 0], "end": [4.05, 2.0], "gpts": 9}),
    ("scanB", abtem.LineScan((0.0, 1.0), (4.05, 1.0), gpts=8, endpoint=False),
     {"start": [0, 1], "end": [4.05, 1.0], "gpts": 8, "endpoint": False}),
]:
    m = probe.scan(pot, scan, abtem.AnnularDetector(inner=50, outer=150)).compute()
    assert type(m).__name__ == "RealSpaceLineProfiles" and m.array.shape == (scan.gpts,)
    put(f"scan.{label}.adf", np.asarray(m.array))
    put(f"scan.{label}.positions", scan.get_positions())
    sampling = float(m.sampling)
    assert sampling == float(scan.sampling)
    entry = {"label": label, **spec, "sampling": sampling}
    # width and interpolate on the profile abTEM returned
    for h in (0.5, 0.2):
        entry[f"width{h}"] = float(np.asarray(m.width(h)))
    for tag, kw in [("g25", dict(gpts=25)), ("g25e", dict(gpts=25, endpoint=True)), ("g31o1", dict(gpts=31, order=1)),
                    ("g31o5", dict(gpts=31, order=5, endpoint=True)), ("s0.1", dict(sampling=0.1)),
                    ("same", dict()), ("g40o0", dict(gpts=40, order=0)), ("g17o2", dict(gpts=17, order=2)),
                    ("g17o4e", dict(gpts=17, order=4, endpoint=True))]:
        r = m.interpolate(**kw)
        put(f"scan.{label}.interp.{tag}", np.asarray(r.array))
        entry.setdefault("interp", {})[tag] = {"args": kw, "sampling": float(r.sampling), "gpts": int(r.array.shape[-1])}
    scan_specs.append(entry)
cases["scans"] = scan_specs

# ---- 3. width and interpolate on planted profiles (ensemble axes, odd and even lengths) -------------------------------
planted = []
for label, shape, sampling in [("p12", (12,), 0.37), ("p3x15", (3, 15), 0.21), ("p2x2x10", (2, 2, 10), 1.3)]:
    x = np.linspace(0, 1, shape[-1], endpoint=False)
    a = np.exp(-((x - 0.45) / 0.18) ** 2) * (1 + rng.random(shape)) + 0.05 * rng.standard_normal(shape)
    p = abtem.RealSpaceLineProfiles(a, sampling=sampling, ensemble_axes_metadata=[UnknownAxis() for _ in shape[:-1]])
    put(f"planted.{label}.array", a)
    entry = {"label": label, "shape": list(shape), "sampling": sampling, "width": {}}
    for h in (0.5, 0.25, 0.8):
        w = np.asarray(p.width(h))
        put(f"planted.{label}.width{h}", np.asarray(w, dtype=float))
        entry["width"][str(h)] = list(w.shape)
    entry["interp"] = {}
    for tag, kw in [("g20", dict(gpts=20)), ("g33e", dict(gpts=33, endpoint=True)), ("g9o5", dict(gpts=9, order=5)),
                    ("s0.05", dict(sampling=0.05)), ("g50o1", dict(gpts=50, order=1)), ("g12o0", dict(gpts=12, order=0)),
                    ("g27o2e", dict(gpts=27, order=2, endpoint=True)), ("g27o4", dict(gpts=27, order=4))]:
        r = p.interpolate(**kw)
        put(f"planted.{label}.interp.{tag}", np.asarray(r.array))
        entry["interp"][tag] = {"args": kw, "sampling": float(r.sampling), "gpts": int(r.array.shape[-1]), "shape": list(r.array.shape)}
    planted.append(entry)
cases["planted"] = planted

# ---- 4. interpolate_line on an Images and on a DiffractionPatterns --------------------------------------------------
img_shape, img_sampling = (24, 18), (0.2, 0.25)
img = rng.standard_normal(img_shape) + 3.0 * np.add.outer(np.sin(np.arange(img_shape[0]) / 3), np.cos(np.arange(img_shape[1]) / 2))
images = abtem.Images(img, sampling=img_sampling)
put("img.array", img)
dp_shape, dp_sampling = (21, 20), (0.11, 0.13)
dpa = rng.random(dp_shape) ** 2
pattern = abtem.DiffractionPatterns(dpa, sampling=dp_sampling, metadata={"energy": 200e3})
put("dp.array", dpa)
cases["image"] = {"shape": list(img_shape), "sampling": list(img_sampling)}
cases["pattern"] = {"shape": list(dp_shape), "sampling": list(dp_sampling), "offset": list(pattern.offset)}

line_specs = []


def record_line(label, obj, kind, **kw):
    r = obj.interpolate_line(**kw) if kind == "line" else obj.interpolate_line_at_position(**kw)
    put(f"{label}.values", np.asarray(r.array))
    line_specs.append({"label": label, "target": "image" if obj is images else "pattern", "kind": kind, "args": kw,
                       "sampling": float(r.sampling), "gpts": int(r.array.shape[-1]),
                       "start": list(r.metadata["start"]), "end": list(r.metadata["end"]),
                       "type": type(r).__name__, "width": r.metadata.get("width")})


for obj, tag in [(images, "img"), (pattern, "dp")]:
    ext = obj.extent
    record_line(f"{tag}.l1", obj, "line", start=(0.5, 0.4), end=(3.9, 3.1), gpts=40)
    record_line(f"{tag}.l2", obj, "line", start=(0.5, 0.4), end=(3.9, 3.1), gpts=40, endpoint=True)
    record_line(f"{tag}.l3", obj, "line", start=(0.5, 0.4), end=(3.9, 3.1), sampling=0.07)
    record_line(f"{tag}.l4", obj, "line", start=(0.5, 0.4), end=(3.9, 3.1), gpts=33, order=1)
    record_line(f"{tag}.l5", obj, "line", start=(0.5, 0.4), end=(3.9, 3.1), gpts=33, order=5)
    record_line(f"{tag}.l6", obj, "line", start=(0.5, 0.4), end=(3.9, 3.1), gpts=33, order=0)
    record_line(f"{tag}.l7", obj, "line", start=(0.5, 0.4), end=(3.9, 3.1), gpts=33, order=2)
    record_line(f"{tag}.l8", obj, "line", start=(0.5, 0.4), end=(3.9, 3.1), gpts=33, order=4)
    record_line(f"{tag}.m1", obj, "line", start=(0.5, 0.4), end=(3.0, 2.5), gpts=30, margin=0.3)
    record_line(f"{tag}.m2", obj, "line", start=(0.5, 0.4), end=(3.0, 2.5), sampling=0.05, margin=0.4, endpoint=True)
    record_line(f"{tag}.w1", obj, "line", start=(0.3, 0.3), end=(3.5, 2.8), gpts=36, width=0.6)
    record_line(f"{tag}.w2", obj, "line", start=(0.3, 0.3), end=(3.5, 2.8), gpts=36, width=0.17, order=1)
    record_line(f"{tag}.w3", obj, "line", start=(0.3, 0.3), end=(3.5, 2.8), gpts=36, width=1.0, margin=0.2, endpoint=True)
    record_line(f"{tag}.f1", obj, "line", start=(0.1, 0.2), end=(0.8, 0.6), gpts=25, fractional=True)
    record_line(f"{tag}.d1", obj, "line", end=(2.5, 3.0), gpts=20)
    record_line(f"{tag}.d2", obj, "line", start=(0.0, 0.0), end=(1.0, 1.0), fractional=True)
    record_line(f"{tag}.a1", obj, "pos", center=(2.0, 1.5), angle=25.0, extent=2.5, gpts=30)
    record_line(f"{tag}.a2", obj, "pos", center=(2.0, 1.5), angle=135.0, extent=2.5, sampling=0.06, width=0.4)
    record_line(f"{tag}.a3", obj, "pos", center=(2.0, 1.5), angle=0.0, extent=2.0, gpts=21, endpoint=False, order=5)
cases["interpLines"] = line_specs

cases["arrays"] = {"file": "linescan-arrays.bin", "sha256": hashlib.sha256(bytes(BLOB)).hexdigest(), "index": INDEX}
with open(os.path.join(OUT, "linescan-arrays.bin"), "wb") as f:
    f.write(BLOB)
with open(os.path.join(OUT, "linescan-pins.json"), "w") as f:
    json.dump(cases, f)
print("ok", len(INDEX), "arrays", len(BLOB), "bytes")
