"""SimulatoKit potential and multislice pins: the pin generator.

Runs abTEM 1.0.10 in float64 (fft numpy, eager) on the cases below, re-derives every quantity from the
published formulas in plain NumPy (the functions named `text_*` below are written from the formulas, not from abTEM's
code paths) and asserts the two agree; then writes

  Tests/SimulatoMultisliceTests/Fixtures/slice-1-pins.json   (scalars, shapes, the array index, versions)
  Tests/SimulatoMultisliceTests/Fixtures/slice-1-arrays.bin  (float64 little-endian; complex as interleaved re, im)
  Sources/SimulatoMultislice/LobatoTable.swift               (abTEM's lobato.json, verbatim numbers, with its sha256)

Run from any directory:  ../References/python/abtem/bin/python -I -B tools/slice-1/pins.py
A failed assertion stops the generator: a `text_*` function then disagrees with abTEM, and the fixtures are not written.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FIXTURES = os.path.join(ROOT, "Tests", "SimulatoMultisliceTests", "Fixtures")
TABLE_SWIFT = os.path.join(ROOT, "Sources", "SimulatoMultislice", "LobatoTable.swift")

import abtem  # noqa: E402
import ase  # noqa: E402
import scipy  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.build import bulk  # noqa: E402
from scipy.special import kn  # noqa: E402

abtem.config.set({"precision": "float64", "fft": "numpy", "dask.lazy": False,
                  "diagnostics.progress_bar": False})

from abtem.core.constants import kappa as ABTEM_KAPPA  # noqa: E402
from abtem.core.energy import energy2sigma, energy2wavelength  # noqa: E402
from abtem.parametrizations import LobatoParametrization  # noqa: E402

LOBATO = LobatoParametrization()
PI32 = float(np.float32(np.pi))
PI2_32 = float(np.float32(np.pi ** 2))
SYMBOLS = ase.data.chemical_symbols

# --------------------------------------------------------------------------------------------------------------------
# The formulas, written independently of abTEM's code.
# --------------------------------------------------------------------------------------------------------------------

H, ME, C, E = 6.62607004e-34, 9.10938356e-31, 299792458.0, 1.6021766208e-19


def text_wavelength(energy):  # 3.2
    return H * C / np.sqrt(energy * (2 * ME * C ** 2 / E + energy)) / E * 1e10


def text_scaled(symbol):  # 3.3
    p = np.array(LOBATO.parameters[symbol], dtype=np.float64)
    a = np.pi ** 2 * p[0] / p[1] ** 1.5 / ABTEM_KAPPA
    b = 2 * np.pi / np.sqrt(p[1])
    return a, b


def text_f(symbol, k2, pi_choice):  # 3.3
    a, b = text_scaled(symbol)
    pi, pi2 = (PI32, PI2_32) if pi_choice == "abTEMFloat32" else (np.pi, np.pi ** 2)
    q = 4 * pi2 * np.asarray(k2, dtype=np.float64)
    out = np.zeros_like(q)
    for i in range(5):  # the order of summation abTEM uses: term by term, i = 0 ... 4
        out = out + (a[i] / b[i] / (q + b[i] ** 2) + a[i] * b[i] / (q + b[i] ** 2) ** 2)
    return 8 * pi * out


def fftfreq2(gpts, sampling):
    kx = np.fft.fftfreq(gpts[0], sampling[0])
    ky = np.fft.fftfreq(gpts[1], sampling[1])
    return kx[:, None], ky[None, :]


def text_compensation(gpts, sampling):  # 3.5 step 3
    kx, ky = fftfreq2(gpts, sampling)
    kappa_ = np.sqrt((kx * sampling[0]) ** 2 + (ky * sampling[1]) ** 2)
    k = kappa_.copy()
    k[0, 0] = 1.0
    s = np.sin(k) / k * sampling[0] * sampling[1]
    s[0, 0] = sampling[0] * sampling[1]
    return s


def text_deltas(xy, gpts, sampling):  # 3.5 step 1
    arr = np.zeros(gpts)
    for x, y in xy:
        u, v = x / sampling[0], y / sampling[1]
        i, j = int(np.floor(u)), int(np.floor(v))
        fx, fy = u - i, v - j
        for di, dj, w in ((0, 0, (1 - fx) * (1 - fy)), (1, 0, fx * (1 - fy)), (0, 1, (1 - fx) * fy), (1, 1, fx * fy)):
            arr[(i + di) % gpts[0], (j + dj) % gpts[1]] += w
    return arr


def text_wrap_and_slice(numbers, positions, box, slice_thickness):  # 3.4
    pos = np.array(positions, dtype=np.float64).reshape(-1, 3).copy()
    for d in range(3):
        pos[:, d] = (pos[:, d] / box[d] % 1.0) * box[d]
    pos[pos[:, 2] > box[2] - 1e-10, 2] = 0.0
    n = int(np.ceil(box[2] / slice_thickness))
    thick = (box[2] / n,) * n
    edges = np.cumsum(thick) - 1e-12
    labels = np.array([np.searchsorted(edges, z, side="right") for z in pos[:, 2]], dtype=int)
    labels = np.minimum(labels, n - 1)
    return pos, thick, labels


def text_potential_slice(numbers, pos, gpts, sampling, pi_choice):  # 3.5
    kx, ky = fftfreq2(gpts, sampling)
    k2 = kx ** 2 + ky ** 2
    s = text_compensation(gpts, sampling)
    v = np.zeros(gpts)
    for z in sorted(set(numbers)):
        xy = [p[:2] for p, zz in zip(pos, numbers) if zz == z]
        d = np.fft.fft2(text_deltas(xy, gpts, sampling).astype(np.complex128))
        v = v + np.fft.ifft2(d * text_f(SYMBOLS[z], k2, pi_choice) / s).real
    return v


def text_antialias(gpts, sampling):  # 3.6
    kc = 0.6666666 / max(sampling) / 2
    w = 0.01 / max(sampling)
    kx, ky = fftfreq2(gpts, sampling)
    r = np.sqrt(kx ** 2 + ky ** 2)
    a = 0.5 * (1 + np.cos(np.pi * (r - kc + w) / w))
    a[r > kc] = 0.0
    return np.where(r > kc - w, a, 1.0)


def text_transmission(v, energy, gpts, sampling):  # 3.6
    t = np.exp(1j * energy2sigma(energy) * v)
    return np.fft.ifft2(np.fft.fft2(t) * text_antialias(gpts, sampling))


def text_propagator(dz, energy, gpts, sampling):  # 3.6
    kx, ky = fftfreq2(gpts, sampling)
    lam = text_wavelength(energy)
    return np.exp(-1j * np.pi * lam * dz * (kx ** 2 + ky ** 2)) * text_antialias(gpts, sampling)


def text_chi(alpha, phi, c, lam):  # 3.8
    g = lambda k: c.get(k, 0.0)  # noqa: E731
    chi = (0.5 * alpha ** 2 * (g("C10") + g("C12") * np.cos(2 * (phi - g("phi12"))))
           + 1 / 3 * alpha ** 3 * (g("C21") * np.cos(phi - g("phi21")) + g("C23") * np.cos(3 * (phi - g("phi23"))))
           + 0.25 * alpha ** 4 * (g("C30") + g("C32") * np.cos(2 * (phi - g("phi32")))
                                  + g("C34") * np.cos(4 * (phi - g("phi34")))))
    return 2 * np.pi / lam * chi


def text_probe(energy, gpts, sampling, cutoff_mrad, aberr, position):  # 3.8
    lam = text_wavelength(energy)
    kx, ky = fftfreq2(gpts, sampling)
    alpha = np.sqrt(kx ** 2 + ky ** 2) * lam
    phi = np.arctan2(ky, kx)
    dax = lam / (gpts[0] * sampling[0])
    day = lam / (gpts[1] * sampling[1])
    den = np.sqrt((np.cos(phi) * dax) ** 2 + (np.sin(phi) * day) ** 2)
    den[0, 0] = 1.0
    ap = np.clip((cutoff_mrad * 1e-3 - alpha) / den + 0.5, 0.0, 1.0)
    ap[0, 0] = 1.0
    ux, uy = position[0] / sampling[0], position[1] / sampling[1]
    kpx = np.fft.fftfreq(gpts[0], 1.0)[:, None]
    kpy = np.fft.fftfreq(gpts[1], 1.0)[None, :]
    shift = np.exp(-2j * np.pi * kpx * ux) * np.exp(-2j * np.pi * kpy * uy)
    psi_k = shift * ap * np.exp(-1j * text_chi(alpha, phi, aberr, lam))
    psi_k = psi_k / np.sqrt((np.abs(psi_k) ** 2).sum())
    return np.fft.ifft2(psi_k)


def text_multislice(psi, potentials, thicknesses, energy, gpts, sampling, exit_after):  # 3.7
    out = []
    for j, (v, dz) in enumerate(zip(potentials, thicknesses)):
        t = text_transmission(v, energy, gpts, sampling)
        psi = np.fft.ifft2(np.fft.fft2(psi * t) * text_propagator(dz, energy, gpts, sampling))
        if j in exit_after:
            out.append(psi)
    return out


def cutoff_gpts_odd(gpts, sampling):  # parity "odd", the diffraction default
    kcut = 2.0 / 3.0 / max(sampling)
    out = []
    for n, d in zip(gpts, sampling):
        m = int(np.floor(np.round(kcut * n * d, 7)))
        if m % 2 != n % 2:  # first to the grid's parity (_antialias_cutoff_gpts, parity "same")
            m += 1
        if m % 2 == 0:  # then odd (diffraction_patterns' default parity "odd")
            m += 1
        out.append(m)
    return tuple(out)


def text_pattern(psi, gpts, sampling, plane_wave):  # 3.10
    a = psi / float(gpts[0] * gpts[1]) if plane_wave else psi
    f = np.fft.fft2(a)
    ng = cutoff_gpts_odd(gpts, sampling)
    keep = []
    for n, m in zip(gpts, ng):
        lo = m // 2 + 1 if m % 2 else m // 2
        hi = m // 2 if m % 2 else m // 2
        keep.append(np.r_[np.arange(lo), np.arange(n - hi, n)])
    f = f[np.ix_(keep[0], keep[1])]
    return np.fft.fftshift(np.abs(f) ** 2, axes=(-1, -2))


# --------------------------------------------------------------------------------------------------------------------
# Output plumbing
# --------------------------------------------------------------------------------------------------------------------

BLOB = bytearray()
INDEX: dict[str, dict] = {}


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


def relmax(a, b):
    a, b = np.asarray(a), np.asarray(b)
    m = np.max(np.abs(b))
    return float(np.max(np.abs(a - b)) / (m if m > 0 else 1.0))


CHECKS: list[str] = []


def check(label, value, bound):
    ok = value <= bound
    CHECKS.append(f"{'ok ' if ok else 'FAIL'} {label}: {value:.3e} <= {bound:.0e}")
    if not ok:
        print("\n".join(CHECKS))
        raise SystemExit(f"self-check failed: {label}")


# --------------------------------------------------------------------------------------------------------------------
# The cases
# --------------------------------------------------------------------------------------------------------------------


def atoms_of(numbers, positions, box):
    return Atoms(numbers=numbers, positions=positions, cell=box, pbc=True)


def al_cubic(nx, ny, nz):
    a = bulk("Al", "fcc", a=4.05, cubic=True) * (nx, ny, nz)
    return a


def srtio3(nx, ny, nz):
    a0 = 3.905
    base = Atoms("SrTiO3", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)],
                 cell=[a0, a0, a0], pbc=True)
    return base * (nx, ny, nz)


def run_case(name, atoms, gpts, slice_thickness, energy, exit_every, beams, store, pattern_store=None, store_potential=True):
    """beams: list of (label, kind, params). store: exit-plane indices whose exit wave is stored in full;
    pattern_store: those whose diffraction pattern is stored (default: store); store_potential: V and t of slices 0, 1."""
    pattern_store = store if pattern_store is None else pattern_store
    box = tuple(float(x) for x in np.diag(atoms.cell))
    sampling = (box[0] / gpts[0], box[1] / gpts[1])
    case = {"box": box, "gpts": list(gpts), "sampling": list(sampling), "energy": energy,
            "slice_thickness_requested": slice_thickness,
            "numbers": [int(z) for z in atoms.numbers], "positions": atoms.positions.tolist()}
    pot = abtem.Potential(atoms, gpts=gpts, slice_thickness=slice_thickness,
                          exit_planes=exit_every if exit_every else None)
    case["slice_thickness"] = list(pot.slice_thickness)
    case["exit_planes"] = [int(e) for e in pot.exit_planes]
    built = pot.build().compute()
    varr = np.asarray(built.array)
    # text re-derivation of the slicing and of every slice (abTEM pi and exact pi)
    pos, thick, labels = text_wrap_and_slice(atoms.numbers, atoms.positions, box, slice_thickness)
    check(f"{name} slice thicknesses", relmax(thick, pot.slice_thickness), 0)
    case["slice_of_atom"] = [int(x) for x in labels]
    v_text_exact = []
    for j in range(len(thick)):
        nums = [int(z) for z, l in zip(atoms.numbers, labels) if l == j]
        pj = [p for p, l in zip(pos, labels) if l == j]
        v32 = text_potential_slice(nums, pj, gpts, sampling, "abTEMFloat32")
        vex = text_potential_slice(nums, pj, gpts, sampling, "exact")
        check(f"{name} V slice {j} text(float32 pi) vs abTEM", relmax(v32, varr[j]), 1e-12)
        check(f"{name} V slice {j} exact pi vs abTEM", relmax(vex, varr[j]), 1e-7)
        if nums:
            area = sampling[0] * sampling[1]
            fsum = sum(text_f(SYMBOLS[z], np.array([0.0]), "abTEMFloat32")[0] for z in nums)
            check(f"{name} V slice {j} k=0 identity", abs(varr[j].sum() * area - fsum) / fsum, 1e-12)
        v_text_exact.append(vex)
    nstore_v = min(2, len(thick)) if store_potential else 0
    for j in range(nstore_v):
        put(f"{name}.V.{j}", varr[j])
    case["V_stored_slices"] = nstore_v
    case["V_max_abs"] = [float(np.max(np.abs(varr[j]))) for j in range(len(thick))]
    case["V_sum"] = [float(varr[j].sum()) for j in range(len(thick))]
    t0 = text_transmission(varr[0], energy, gpts, sampling)
    from abtem.antialias import AntialiasAperture
    tf = built.transmission_function(energy)
    t0ab = np.asarray(AntialiasAperture().bandlimit(tf, in_place=False).array)[0]
    check(f"{name} bandlimited t[0] text vs abTEM", relmax(t0, t0ab), 1e-12)
    if store_potential:
        put(f"{name}.t.0", t0)

    exit_after = [e for e in pot.exit_planes if e >= 0]
    case["beams"] = {}
    for label, kind, params in beams:
        b = {"kind": kind, **params}
        if kind == "plane":
            wave = abtem.PlaneWave(energy=energy)
            incident = np.ones(gpts, dtype=np.complex128)
            ex = wave.multislice(pot).compute()
        else:
            aberr = params["aberrations"]
            probe = abtem.Probe(energy=energy, semiangle_cutoff=params["cutoff_mrad"], gpts=gpts, extent=box[:2],
                                **aberr)
            built_probe = probe.build(scan=[tuple(params["position"])]).compute()
            incident = np.asarray(built_probe.array).reshape(gpts)
            text_inc = text_probe(energy, gpts, sampling, params["cutoff_mrad"], aberr, params["position"])
            check(f"{name}/{label} incident probe text vs abTEM", relmax(text_inc, incident), 1e-12)
            b["incident_sum_k"] = float((np.abs(np.fft.fft2(incident)) ** 2).sum())
            ex = probe.multislice(pot, scan=[tuple(params["position"])]).compute()
        exits = np.asarray(ex.array).reshape((-1,) + tuple(gpts))
        planes = case["exit_planes"]
        assert exits.shape[0] == len(planes), (exits.shape, planes)
        # text multislice from the abTEM potential (so the propagation is checked on its own)
        text_exits = text_multislice(incident, list(varr), pot.slice_thickness, energy, gpts, sampling, exit_after)
        offset = 1 if planes[0] == -1 else 0
        for q, te in enumerate(text_exits):
            check(f"{name}/{label} exit plane {planes[q + offset]} text vs abTEM", relmax(te, exits[q + offset]),
                  1e-11)
        pats = ex.diffraction_patterns(max_angle="cutoff").compute()
        parr = np.asarray(pats.array).reshape((exits.shape[0],) + pats.array.shape[-2:])
        ints = np.abs(exits) ** 2
        for q in range(exits.shape[0]):
            tp = text_pattern(exits[q], gpts, sampling, kind == "plane")
            check(f"{name}/{label} pattern {q} text vs abTEM", relmax(tp, parr[q]), 1e-11)
        b["pattern_gpts"] = list(parr.shape[-2:])
        b["total_intensity"] = [float(ints[q].sum()) for q in range(exits.shape[0])]
        b["pattern_sum"] = [float(parr[q].sum()) for q in range(exits.shape[0])]
        b["exit_max_abs"] = [float(np.abs(exits[q]).max()) for q in range(exits.shape[0])]
        c0, c1 = parr.shape[-2] // 2, parr.shape[-1] // 2
        b["pattern_center"] = [float(parr[q][c0, c1]) for q in range(exits.shape[0])]
        stored = sorted(set(s if s >= 0 else exits.shape[0] + s for s in store))
        pstored = sorted(set(s if s >= 0 else exits.shape[0] + s for s in pattern_store))
        b["stored_planes"] = stored
        b["stored_pattern_planes"] = pstored
        if kind == "probe":
            put(f"{name}/{label}.incident", incident)
        for q in stored:
            put(f"{name}/{label}.exit.{q}", exits[q])
        for q in pstored:
            put(f"{name}/{label}.pattern.{q}", parr[q])
        case["beams"][label] = b
    return case


def main():
    versions = {"abtem": abtem.__version__, "ase": ase.__version__, "numpy": np.__version__,
                "scipy": scipy.__version__, "python": platform.python_version(),
                "abtem_commit": "164e644ff411167db23ec4ef26aeb4fe8f384b26",
                "config": {"precision": "float64", "fft": "numpy", "dask.lazy": False}}
    pins: dict = {"versions": versions, "kappa": ABTEM_KAPPA, "pi_float32": PI32, "pi2_float32": PI2_32}

    # 3.2: wavelength and sigma
    energies = [60e3, 80e3, 100e3, 120e3, 200e3, 300e3]
    pins["beam"] = [{"energy": e, "wavelength": energy2wavelength(e), "sigma": energy2sigma(e)} for e in energies]
    for row in pins["beam"]:
        check(f"lambda {row['energy']:.0f} text vs abTEM", abs(text_wavelength(row["energy"]) - row["wavelength"])
              / row["wavelength"], 1e-15)

    # 3.3: f(k2) for a set of elements, both pi choices
    k2 = np.array([0.0, 0.01, 0.1, 0.25, 1.0, 4.0, 16.0, 100.0])
    pins["f_k2"] = k2.tolist()
    pins["f"] = {}
    for sym in ("H", "C", "O", "Al", "Si", "Ti", "Sr", "Au", "U"):
        fab = LOBATO.projected_scattering_factor(sym)(k2)
        f32 = text_f(sym, k2, "abTEMFloat32")
        fex = text_f(sym, k2, "exact")
        check(f"f {sym} text(float32 pi) vs abTEM", float(np.max(np.abs(f32 - fab) / np.abs(fab))), 1e-14)
        check(f"f {sym} exact pi vs abTEM", float(np.max(np.abs(fex - fab) / np.abs(fab))), 1e-7)
        pins["f"][sym] = {"abtem": fab.tolist(), "exact": fex.tolist()}

    cases = {}
    cases["A1"] = run_case("A1", atoms_of([79], [(5.0, 5.0, 1.0)], (10.0, 10.0, 2.0)), (128, 128), 2.0, 200e3, None,
                           [("plane", "plane", {})], store=[0])
    cases["A2"] = run_case("A2", atoms_of([6, 14], [(4.3, 6.1, 0.5), (12.7, 3.2, 1.5)], (16.0, 16.0, 2.0)),
                           (160, 160), 2.0, 300e3, None, [("plane", "plane", {})], store=[0])
    b2_aberr = {"C10": -50.0, "C30": 1.0e4, "C12": 20.0, "phi12": 0.3}
    al = al_cubic(4, 4, 10)
    cases["B1"] = run_case("B1", al, (256, 256), 2.025, 200e3, 2, [("plane", "plane", {})], store=[-1],
                           pattern_store=[1, 5, -1])
    cases["B2"] = run_case("B2", al, (256, 256), 2.025, 200e3, 2,
                           [("probe", "probe", {"cutoff_mrad": 20.0, "aberrations": b2_aberr,
                                                 "position": [8.37, 7.91]})], store=[-1],
                           store_potential=False)
    vac = Atoms(cell=(16.2, 16.2, 20.0), pbc=True)
    cases["C1"] = run_case("C1", vac, (256, 256), 2.0, 200e3, 1,
                           [("plane", "plane", {}),
                            ("probe", "probe", {"cutoff_mrad": 20.0, "aberrations": b2_aberr, "position": [8.1, 8.1]})],
                           store=[-1], store_potential=False)
    sto = srtio3(3, 4, 5)
    cases["D1"] = run_case("D1", sto, (120, 160), 1.9525, 300e3, None,
                           [("plane", "plane", {}),
                            ("probe", "probe", {"cutoff_mrad": 25.0, "aberrations": {}, "position": [5.3, 7.9]})],
                           store=[0])
    pins["cases"] = cases

    # Truth 3: the analytic projected potential of one Al atom (SciPy kn) against the FFT method, on a
    # grid point at the centre of a 25.6 A box, 512^2 (0.05 A) and 1024^2 (0.025 A), r = 0.25 ... 4 A along x.
    a_, b_ = text_scaled("Al")
    ff = {"box": [25.6, 25.6, 2.0], "position": [12.8, 12.8, 1.0], "grids": {}}
    errs = {}
    for n in (512, 1024):
        d = 25.6 / n
        offs = np.arange(int(round(0.25 / d)), int(round(4.0 / d)) + 1)
        r = offs * d
        vr = 2 * sum(2 * a_[i] / b_[i] * kn(0, r * b_[i]) + a_[i] * r * kn(1, r * b_[i]) for i in range(5))
        vfft = text_potential_slice([13], [np.array([12.8, 12.8, 1.0])], (n, n), (d, d), "exact")
        line = vfft[n // 2 + offs, n // 2]
        near = r <= 2.0 + 1e-9
        rel_near = float(np.max(np.abs(line[near] - vr[near]) / vr[near]))
        abs_far = float(np.max(np.abs(line[~near] - vr[~near])))
        errs[n] = float(np.max(np.abs(line - vr)))
        ff["grids"][str(n)] = {"pixel_offsets": offs.tolist(), "r": r.tolist(), "V": vr.tolist(),
                              "generator_rel_near": rel_near, "generator_abs_far": abs_far,
                              "generator_max_abs": errs[n]}
        if n == 512:
            check("truth 3a far field 512, rel for r <= 2 A (generator)", rel_near, 2e-2)
            check("truth 3b far field 512, abs for 2 < r <= 4 A (generator)", abs_far, 3e-2)
    check("truth 3c halving the pixel halves the max abs error (generator)", errs[1024] / errs[512], 0.5)
    pins["truth_far_field_Al"] = ff

    os.makedirs(FIXTURES, exist_ok=True)
    blob = bytes(BLOB)
    with open(os.path.join(FIXTURES, "slice-1-arrays.bin"), "wb") as fh:
        fh.write(blob)
    pins["arrays"] = {"file": "slice-1-arrays.bin", "sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob),
                      "index": INDEX}
    pins["self_checks"] = CHECKS
    with open(os.path.join(FIXTURES, "slice-1-pins.json"), "w") as fh:
        json.dump(pins, fh, indent=1, sort_keys=True)

    write_table()
    print("\n".join(CHECKS))
    print(f"wrote {len(INDEX)} arrays, {len(blob) / 1e6:.2f} MB; {len(CHECKS)} self-checks passed")


def write_table():
    path = os.path.join(os.path.dirname(abtem.parametrizations.__file__), "data", "lobato.json")
    raw = open(path, "rb").read()
    data = json.loads(raw)
    lines = [
        "//",
        "//  LobatoTable.swift",
        "//  Role: the Lobato & Van Dyck (2014) parameters (a_i, b_i), i = 1 ... 5, per element Z = 1 ... 103, as abTEM",
        "//        1.0.10 ships them (abtem/parametrizations/data/lobato.json, sha256 " + hashlib.sha256(raw).hexdigest()[:16]
        + "...).",
        "//        Generated by tools/slice-1/pins.py; do not edit. Numbers verbatim (Python repr of each double).",
        "//        I. Lobato, D. Van Dyck, Acta Cryst. A 70, 636-649 (2014). abTEM is GPL-3.0-or-later (NOTICE).",
        "//",
        "",
        "/// `lobatoTable[Z - 1]` = ((a1 ... a5), (b1 ... b5)); `nil` where abTEM's table has no entry.",
        "let lobatoTable: [(a: [Double], b: [Double])?] = [",
    ]
    for z in range(1, 104):
        sym = SYMBOLS[z]
        if sym in data:
            a, b = data[sym]
            aa = ", ".join(repr(float(x)) for x in a)
            bb = ", ".join(repr(float(x)) for x in b)
            lines.append(f"    (a: [{aa}], b: [{bb}]), // {z} {sym}")
        else:
            lines.append(f"    nil, // {z} {sym}")
    lines += ["]", "", f"/// sha256 of abTEM's lobato.json the table was generated from.",
              f"let lobatoTableSourceSHA256 = \"{hashlib.sha256(raw).hexdigest()}\"", ""]
    os.makedirs(os.path.dirname(TABLE_SWIFT), exist_ok=True)
    with open(TABLE_SWIFT, "w") as fh:
        fh.write("\n".join(lines))


if __name__ == "__main__":
    main()
