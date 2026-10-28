"""SimulatoKit frozen phonons, scans and detectors: the pin generator.

Runs abTEM 1.0.10 in float64 (numpy FFT, eager) on these cases: frozen phonons (seeds, displacements, the
potential of each configuration), GridScan positions, and the annular, pixelated and segmented detectors averaged over
configurations. Every quantity is re-derived from the published formulas in NumPy (the `text_*` functions here and in
tools/slice-1/pins.py) and asserted equal to abTEM's; a failure stops the generator. Writes

  Tests/SimulatoDetectorTests/Fixtures/slice-2-pins.json
  Tests/SimulatoDetectorTests/Fixtures/slice-2-arrays.bin   (float64 little-endian; complex interleaved)

Run:  ../References/python/abtem/bin/python -I -B tools/slice-2/pins.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIXTURES = os.path.join(ROOT, "Tests", "SimulatoDetectorTests", "Fixtures")

spec = importlib.util.spec_from_file_location("s1", os.path.join(ROOT, "tools", "slice-1", "pins.py"))
s1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s1)  # sets abTEM's configuration (float64, numpy FFT, eager)

import abtem  # noqa: E402
from abtem.inelastic.phonons import validate_seeds  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.build import bulk  # noqa: E402

SYMBOLS = s1.SYMBOLS
BLOB = bytearray()
INDEX: dict = {}
CHECKS: list[str] = []


def put(name, array):
    a = np.ascontiguousarray(array)
    if np.iscomplexobj(a):
        raw = np.stack([a.real, a.imag], axis=-1).astype("<f8").tobytes()
        kind = "c16"
    else:
        raw = a.astype("<f8").tobytes()
        kind = "f8"
    assert name not in INDEX, name
    INDEX[name] = {"offset": len(BLOB), "count": int(a.size), "shape": list(a.shape), "kind": kind}
    BLOB.extend(raw)


def check(label, value, bound):
    ok = value <= bound
    CHECKS.append(f"{'ok ' if ok else 'FAIL'} {label}: {value:.3e} <= {bound:.0e}")
    if not ok:
        print("\n".join(CHECKS))
        raise SystemExit(f"self-check failed: {label}")


relmax = s1.relmax

# --------------------------------------------------------------------------------------------------------------------
# The formulas, written independently of abTEM's code
# --------------------------------------------------------------------------------------------------------------------


def text_displacements(numbers, sigmas, seed, directions="xyz"):  # 3.1
    rng = np.random.default_rng(seed)
    r = rng.normal(size=(len(numbers), 3))
    sig = np.array([sigmas.get(SYMBOLS[z], 0.0) for z in numbers])
    d = np.zeros((len(numbers), 3))
    for axis, name in enumerate("xyz"):
        if name in directions:
            d[:, axis] = sig * r[:, axis]
    return d


def text_positions(start, end, gpts, endpoint=False):  # 3.3
    xs = [np.linspace(s, e, n, endpoint=endpoint) for s, e, n in zip(start, end, gpts)]
    return np.stack(np.meshgrid(*xs, indexing="ij"), axis=-1)


def angular_sampling(gpts, sampling, lam):  # mrad per pixel
    return tuple(1.0 / (n * d) * lam * 1e3 for n, d in zip(gpts, sampling))


def cutoff_gpts_same(gpts, sampling):  # the potential pins' cutoff_gpts with parity "same"
    kcut = 2.0 / 3.0 / max(sampling)
    out = []
    for n, d in zip(gpts, sampling):
        m = int(np.floor(np.round(kcut * n * d, 7)))
        if m % 2 != n % 2:
            m += 1
        out.append(m)
    return tuple(out)


def crop_shift(f, gpts, new):  # fft_crop to `new`, then fftshift
    keep = []
    for n, m in zip(gpts, new):
        if m >= n:
            keep.append(np.arange(n))
            continue
        lo = m // 2 + 1 if m % 2 else m // 2
        hi = m // 2
        keep.append(np.r_[np.arange(lo), np.arange(n - hi, n)])
    return np.fft.fftshift(f[np.ix_(keep[0], keep[1])], axes=(-1, -2))


def text_annular(psi, gpts, sampling, lam, inner, outer):  # 3.4
    p = np.abs(np.fft.fft2(psi)) ** 2
    da = angular_sampling(gpts, sampling, lam)
    if outer is None:
        cg = cutoff_gpts_same(gpts, sampling)
        outer = np.floor(min(cg[0] // 2 * da[0], cg[1] // 2 * da[1]))
    ax = np.fft.fftfreq(gpts[0], 1 / da[0] / gpts[0])[:, None]
    ay = np.fft.fftfreq(gpts[1], 1 / da[1] / gpts[1])[None, :]
    a2 = ax ** 2 + ay ** 2
    mask = (a2 >= inner ** 2) & (a2 < outer ** 2)
    return float((p * mask).sum())


def text_pixelated(psi, gpts, sampling, mode):  # 3.4
    if mode == "full":
        new = gpts
    else:
        cg = cutoff_gpts_same(gpts, sampling)
        if mode == "cutoff":
            new = cg
        else:  # valid
            new = []
            for n, c in zip(gpts, cg):
                m = int(np.floor(np.round(c / np.sqrt(2), 7)))
                if m % 2 != n % 2:
                    m += 1
                new.append(m)
            new = tuple(new)
    return crop_shift(np.abs(np.fft.fft2(psi)) ** 2, gpts, new)


def text_segmented(psi, gpts, sampling, lam, nr, nphi, inner, outer, rotation):  # 3.4
    da = angular_sampling(gpts, sampling, lam)
    new = []
    for n, d in zip(gpts, da):
        m = int(2 * np.ceil(outer / d)) + 1
        if m % 2 != n % 2:
            m += 1
        new.append(m)
    pat = crop_shift(np.abs(np.fft.fft2(psi)) ** 2, gpts, new)
    ax = np.fft.fftfreq(new[0], 1 / da[0] / new[0])[:, None]
    ay = np.fft.fftfreq(new[1], 1 / da[1] / new[1])[None, :]
    alpha = np.sqrt(ax ** 2 + ay ** 2)
    phi = np.arctan2(ay + 0 * ax, ax + 0 * ay)
    phi = (phi - rotation) % (2 * np.pi)
    valid = (alpha >= inner) & (alpha < outer)
    rb = np.floor(nr * (alpha - inner) / (outer - inner)).astype(int)
    ab = np.clip(np.floor(nphi * (phi / (2 * np.pi))), 0, nphi - 1).astype(int)
    bins = np.where(valid, ab + rb * nphi, -1)
    bins = np.fft.fftshift(bins)
    out = np.zeros(nr * nphi)
    for b in range(nr * nphi):
        out[b] = pat[bins == b].sum()
    return out.reshape(nr, nphi)


def config_potentials(atoms, disp_list, gpts, slice_thickness, pi_choice="abTEMFloat32"):
    box = tuple(float(x) for x in np.diag(atoms.cell))
    sampling = (box[0] / gpts[0], box[1] / gpts[1])
    out = []
    for d in disp_list:
        pos, thick, labels = s1.text_wrap_and_slice(atoms.numbers, atoms.positions + d, box, slice_thickness)
        slices = []
        for j in range(len(thick)):
            nums = [int(z) for z, l in zip(atoms.numbers, labels) if l == j]
            pj = [p for p, l in zip(pos, labels) if l == j]
            slices.append(s1.text_potential_slice(nums, pj, gpts, sampling, pi_choice))
        out.append((slices, thick))
    return out


# --------------------------------------------------------------------------------------------------------------------
# Cases
# --------------------------------------------------------------------------------------------------------------------

DETECTORS = [  # (name, kind, params)
    ("BF", "annular", {"inner": 0.0, "outer": 20.0}),
    ("ADF", "annular", {"inner": 30.0, "outer": 60.0}),
    ("HAADF", "annular", {"inner": 70.0, "outer": 130.0}),
    ("annular_default_outer", "annular", {"inner": 0.0, "outer": None}),
    ("pixelated_cutoff", "pixelated", {"max_angle": "cutoff"}),
    ("pixelated_valid", "pixelated", {"max_angle": "valid"}),
    ("pixelated_full", "pixelated", {"max_angle": "full"}),
    ("DPC", "segmented", {"nbins_radial": 2, "nbins_azimuthal": 4, "inner": 0.0, "outer": 20.0, "rotation": 0.0}),
    ("DPC_rotated", "segmented", {"nbins_radial": 2, "nbins_azimuthal": 4, "inner": 0.0, "outer": 20.0,
                                  "rotation": 0.3}),
]


def abtem_detector(kind, params):
    if kind == "annular":
        return abtem.AnnularDetector(inner=params["inner"], outer=params["outer"])
    if kind == "pixelated":
        return abtem.PixelatedDetector(max_angle=params["max_angle"])
    return abtem.SegmentedDetector(**params)


def stem_case(name, atoms, gpts, slice_thickness, energy, sigmas, num_configs, seed, cutoff, scan, store_positions):
    box = tuple(float(x) for x in np.diag(atoms.cell))
    sampling = (box[0] / gpts[0], box[1] / gpts[1])
    lam = s1.energy2wavelength(energy)
    case = {"box": box, "gpts": list(gpts), "sampling": list(sampling), "energy": energy,
            "slice_thickness_requested": slice_thickness, "numbers": [int(z) for z in atoms.numbers],
            "positions": atoms.positions.tolist(), "sigmas": sigmas, "num_configs": num_configs, "seed": seed,
            "cutoff_mrad": cutoff, "scan": scan}
    if sigmas:
        fp = abtem.FrozenPhonons(atoms, num_configs=num_configs, sigmas=sigmas, seed=seed)
        seeds = [int(s) for s in validate_seeds(seed, num_seeds=num_configs)]
        assert tuple(seeds) == tuple(int(s) for s in fp.seed), (seeds, fp.seed)
        source = fp
    else:
        seeds = []
        source = atoms
    case["seeds"] = seeds
    disp = [text_displacements(atoms.numbers, sigmas, s) for s in seeds] if sigmas else [np.zeros((len(atoms), 3))]
    for c, d in enumerate(disp):
        put(f"{name}.displacement.{c}", d)
    pot = abtem.Potential(source, gpts=gpts, slice_thickness=slice_thickness)
    built = np.asarray(pot.build().compute().array)
    if built.ndim == 3:
        built = built[None]
    tp = config_potentials(atoms, disp, gpts, slice_thickness)
    for c, (slices, thick) in enumerate(tp):
        for j, v in enumerate(slices):
            check(f"{name} config {c} V slice {j} text vs abTEM", relmax(v, built[c, j]), 1e-12)
    put(f"{name}.V.config0.slice0", built[0, 0])
    case["slice_thickness"] = list(tp[0][1])
    case["num_slices"] = len(tp[0][1])

    gs = abtem.GridScan(start=scan["start"], end=scan["end"], gpts=scan["gpts"], endpoint=False)
    positions = np.asarray(gs.get_positions())
    check(f"{name} scan positions text vs abTEM", relmax(text_positions(scan["start"], scan["end"], scan["gpts"]),
                                                          positions), 0)
    put(f"{name}.scan_positions", positions)
    probe = abtem.Probe(energy=energy, semiangle_cutoff=cutoff)
    dets = [abtem_detector(k, p) for _, k, p in DETECTORS]
    measured = probe.scan(pot, scan=gs, detectors=dets)
    measured = [np.asarray(m.compute().array) for m in measured]

    # text: every configuration, every position
    flat = positions.reshape(-1, 2)
    text = {n: [] for n, _, _ in DETECTORS}
    per_config_haadf = []
    for c, (slices, thick) in enumerate(tp):
        vals = {n: [] for n, _, _ in DETECTORS}
        for pos in flat:
            psi0 = s1.text_probe(energy, gpts, sampling, cutoff, {}, pos)
            psi = s1.text_multislice(psi0, list(built[c]), thick, energy, gpts, sampling, [len(thick) - 1])[0]
            for n, k, p in DETECTORS:
                if k == "annular":
                    vals[n].append(text_annular(psi, gpts, sampling, lam, p["inner"], p["outer"]))
                elif k == "pixelated":
                    vals[n].append(text_pixelated(psi, gpts, sampling, p["max_angle"]))
                else:
                    vals[n].append(text_segmented(psi, gpts, sampling, lam, p["nbins_radial"], p["nbins_azimuthal"],
                                                  p["inner"], p["outer"], p["rotation"]))
        per_config_haadf.append(vals["HAADF"][0])
        for n in text:
            text[n].append(np.array(vals[n]))
    case["per_config_HAADF_position0"] = per_config_haadf
    case["detectors"] = {}
    sg = tuple(scan["gpts"])
    for (n, k, p), m in zip(DETECTORS, measured):
        mean = np.mean(np.stack(text[n]), axis=0)  # 3.5: the mean over configurations
        mean = mean.reshape(sg + mean.shape[1:])
        m = m.reshape(sg + mean.shape[len(sg):])
        # abTEM accumulates the segments in float32 (measurements.py 3056-3061, numba fastmath), even in float64 mode:
        # the segmented values agree only to float32 rounding; the text's float64 sum is pinned beside it.
        # Relative to the larger of the image's own maximum and the total intensity (the default-outer annular image's
        # maximum): an image that is pure rounding noise (vacuum ADF beyond the aperture) is judged on the total's scale.
        scale = max(float(np.max(np.abs(m))), float(np.max(np.mean(np.stack(text["annular_default_outer"]), axis=0))))
        diff = float(np.max(np.abs(mean - m))) / scale
        check(f"{name} {n} text vs abTEM", diff, 1e-6 if k == "segmented" else 1e-11)
        entry_scale = scale
        entry = {"kind": k, "params": p, "shape": list(m.shape), "scale": entry_scale}
        if k == "segmented":
            put(f"{name}.{n}.float64", mean)
        if k == "pixelated":
            for q in store_positions:
                put(f"{name}.{n}.position.{q[0]}.{q[1]}", m[q[0], q[1]])
            entry["stored_positions"] = store_positions
            entry["sums"] = m.reshape(sg + (-1,)).sum(-1).tolist()
        else:
            put(f"{name}.{n}", m)
        case["detectors"][n] = entry
    return case


def frozen_plane_case(name, atoms, gpts, slice_thickness, energy, sigmas, num_configs, seed):
    box = tuple(float(x) for x in np.diag(atoms.cell))
    sampling = (box[0] / gpts[0], box[1] / gpts[1])
    fp = abtem.FrozenPhonons(atoms, num_configs=num_configs, sigmas=sigmas, seed=seed)
    seeds = [int(s) for s in fp.seed]
    disp = [text_displacements(atoms.numbers, sigmas, s) for s in seeds]
    pot = abtem.Potential(fp, gpts=gpts, slice_thickness=slice_thickness)
    built = np.asarray(pot.build().compute().array)
    exits = np.asarray(abtem.PlaneWave(energy=energy).multislice(pot).compute().array)
    pats = np.asarray(abtem.PlaneWave(energy=energy).multislice(pot).diffraction_patterns(max_angle="cutoff")
                      .compute().array)
    tp = config_potentials(atoms, disp, gpts, slice_thickness)
    texts = []
    for c, (slices, thick) in enumerate(tp):
        for j, v in enumerate(slices):
            check(f"{name} config {c} V slice {j} text vs abTEM", relmax(v, built[c, j]), 1e-12)
        psi = s1.text_multislice(np.ones(gpts, complex), list(built[c]), thick, energy, gpts, sampling,
                                 [len(thick) - 1])[0]
        check(f"{name} config {c} exit text vs abTEM", relmax(psi, exits[c]), 1e-11)
        texts.append(s1.text_pattern(psi, gpts, sampling, True))
    mean = np.mean(np.stack(texts), axis=0)
    check(f"{name} averaged pattern text vs abTEM", relmax(mean, pats.mean(0)), 1e-11)
    for c, d in enumerate(disp):
        put(f"{name}.displacement.{c}", d)
    put(f"{name}.V.config0.slice0", built[0, 0])
    put(f"{name}.pattern_mean", pats.mean(0))
    return {"box": box, "gpts": list(gpts), "energy": energy, "slice_thickness_requested": slice_thickness,
            "numbers": [int(z) for z in atoms.numbers], "positions": atoms.positions.tolist(), "sigmas": sigmas,
            "num_configs": num_configs, "seed": seed, "seeds": seeds, "slice_thickness": list(tp[0][1]),
            "pattern_shape": list(pats.shape[1:]), "per_config_total": [float((np.abs(e) ** 2).sum()) for e in exits]}


def main():
    pins = {"versions": {"abtem": abtem.__version__, "numpy": np.__version__, "python": sys.version.split()[0],
                         "abtem_commit": "164e644ff411167db23ec4ef26aeb4fe8f384b26"},
            "detectors": [[n, k, p] for n, k, p in DETECTORS], "cases": {}}
    a0 = 3.905
    sto = Atoms("SrTiO3", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)],
                cell=[a0, a0, a0], pbc=True) * (2, 2, 4)
    pins["cases"]["F1"] = frozen_plane_case("F1", sto, (120, 120), 1.9525, 200e3,
                                            {"Sr": 0.08, "Ti": 0.06, "O": 0.07}, 4, 7)
    al = bulk("Al", "fcc", a=4.05, cubic=True) * (2, 2, 3)
    scan = {"start": [0.0, 0.0], "end": [4.05, 4.05], "gpts": [6, 6]}
    pins["cases"]["S1"] = stem_case("S1", al, (128, 128), 2.025, 200e3, {"Al": 0.085}, 3, 3, 20.0, scan,
                                    [[0, 0], [2, 3]])
    pins["cases"]["S2"] = stem_case("S2", al, (128, 128), 2.025, 200e3, {}, 1, None, 20.0, scan, [[0, 0]])
    vac = Atoms(cell=(8.1, 8.1, 6.0), pbc=True)
    pins["cases"]["S3"] = stem_case("S3", vac, (128, 128), 2.0, 200e3, {}, 1, None, 20.0, scan, [[0, 0]])
    single = {}
    for sym in ("C", "Si", "Cu", "Au"):
        atoms = Atoms(sym, positions=[(4.0, 4.0, 1.0)], cell=(8.0, 8.0, 2.0), pbc=True)
        pot = abtem.Potential(atoms, gpts=(128, 128), slice_thickness=2.0)
        probe = abtem.Probe(energy=200e3, semiangle_cutoff=25.0)
        ex = probe.multislice(pot, scan=[(4.0, 4.0)]).compute()
        haadf = float(np.asarray(abtem.AnnularDetector(inner=70, outer=130).detect(ex).compute().array).ravel()[0])
        psi = np.asarray(ex.array).reshape(128, 128)
        t = text_annular(psi, (128, 128), (8 / 128, 8 / 128), s1.energy2wavelength(200e3), 70, 130)
        check(f"S4 {sym} HAADF text vs abTEM", abs(t - haadf) / haadf, 1e-11)
        single[sym] = haadf
    pins["cases"]["S4"] = {"HAADF": single, "box": [8.0, 8.0, 2.0], "gpts": [128, 128], "energy": 200e3,
                           "cutoff_mrad": 25.0, "position": [4.0, 4.0, 1.0]}

    os.makedirs(FIXTURES, exist_ok=True)
    blob = bytes(BLOB)
    with open(os.path.join(FIXTURES, "slice-2-arrays.bin"), "wb") as fh:
        fh.write(blob)
    pins["arrays"] = {"file": "slice-2-arrays.bin", "sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob),
                      "index": INDEX}
    pins["self_checks"] = CHECKS
    with open(os.path.join(FIXTURES, "slice-2-pins.json"), "w") as fh:
        json.dump(pins, fh, indent=1, sort_keys=True)
    print("\n".join(CHECKS))
    print(f"wrote {len(INDEX)} arrays, {len(blob) / 1e6:.2f} MB; {len(CHECKS)} self-checks passed")


if __name__ == "__main__":
    main()
