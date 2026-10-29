"""SimulatoKit CTF and HRTEM images: the pin generator.

abTEM 1.0.10 (float64, numpy FFT, eager): CTF kernels, HRTEM images (apply_ctf then intensity) of a SrTiO3 exit wave, and
the closed forms. Every quantity is re-derived from the published formulas in NumPy and asserted equal to abTEM's.
Writes Tests/SimulatoKitParityTests/Fixtures/slice-3-pins.json and slice-3-arrays.bin.

Run from the repo root:  .venv/bin/python -I -B tools/pins/ctf/pins.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FIXTURES = os.path.join(ROOT, "Tests", "SimulatoKitParityTests", "Fixtures")
spec = importlib.util.spec_from_file_location("s1", os.path.join(ROOT, "tools", "pins", "potential-multislice", "pins.py"))
s1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s1)

import abtem  # noqa: E402
from abtem.transfer import point_resolution, scherzer_defocus, nyquist_sampling  # noqa: E402
from ase import Atoms  # noqa: E402

BLOB = bytearray()
INDEX: dict = {}
CHECKS: list[str] = []


def put(name, a):
    a = np.ascontiguousarray(a)
    if np.iscomplexobj(a):
        raw, kind = np.stack([a.real, a.imag], axis=-1).astype("<f8").tobytes(), "c16"
    else:
        raw, kind = a.astype("<f8").tobytes(), "f8"
    INDEX[name] = {"offset": len(BLOB), "count": int(a.size), "shape": list(a.shape), "kind": kind}
    BLOB.extend(raw)


def check(label, value, bound):
    ok = value <= bound
    CHECKS.append(f"{'ok ' if ok else 'FAIL'} {label}: {value:.3e} <= {bound:.0e}")
    if not ok:
        print("\n".join(CHECKS))
        raise SystemExit(f"self-check failed: {label}")


SYM = ["C10", "C12", "phi12", "C21", "phi21", "C23", "phi23", "C30", "C32", "phi32", "C34", "phi34", "C41", "phi41",
       "C43", "phi43", "C45", "phi45", "C50", "C52", "phi52", "C54", "phi54", "C56", "phi56"]


def text_chi(alpha, phi, c, lam):  # the aberration phase χ, all orders
    g = lambda k: c.get(k, 0.0)  # noqa: E731
    cs = np.cos
    x = (0.5 * alpha ** 2 * (g("C10") + g("C12") * cs(2 * (phi - g("phi12"))))
         + 1 / 3 * alpha ** 3 * (g("C21") * cs(phi - g("phi21")) + g("C23") * cs(3 * (phi - g("phi23"))))
         + 0.25 * alpha ** 4 * (g("C30") + g("C32") * cs(2 * (phi - g("phi32"))) + g("C34") * cs(4 * (phi - g("phi34"))))
         + 0.2 * alpha ** 5 * (g("C41") * cs(phi - g("phi41")) + g("C43") * cs(3 * (phi - g("phi43")))
                               + g("C45") * cs(5 * (phi - g("phi45"))))
         + 1 / 6 * alpha ** 6 * (g("C50") + g("C52") * cs(2 * (phi - g("phi52"))) + g("C54") * cs(4 * (phi - g("phi54")))
                                 + g("C56") * cs(6 * (phi - g("phi56")))))
    return 2 * np.pi / lam * x


def text_spatial(alpha, phi, c, lam, beta_mrad):  # 3.2
    g = lambda k: c.get(k, 0.0)  # noqa: E731
    cs, sn = np.cos, np.sin
    b = beta_mrad / 1e3
    dk = 2 * np.pi / lam * ((g("C12") * cs(2 * (phi - g("phi12"))) + g("C10")) * alpha
                            + (g("C23") * cs(3 * (phi - g("phi23"))) + g("C21") * cs(phi - g("phi21"))) * alpha ** 2
                            + (g("C34") * cs(4 * (phi - g("phi34"))) + g("C32") * cs(2 * (phi - g("phi32"))) + g("C30"))
                            * alpha ** 3
                            + (g("C45") * cs(5 * (phi - g("phi45"))) + g("C43") * cs(3 * (phi - g("phi43")))
                               + g("C41") * cs(phi - g("phi41"))) * alpha ** 4
                            + (g("C56") * cs(6 * (phi - g("phi56"))) + g("C54") * cs(4 * (phi - g("phi54")))
                               + g("C52") * cs(2 * (phi - g("phi52"))) + g("C50")) * alpha ** 5)
    dp = -2 * np.pi / lam * (0.5 * (2 * g("C12") * sn(2 * (phi - g("phi12")))) * alpha
                             + 1 / 3 * (3 * g("C23") * sn(3 * (phi - g("phi23"))) + g("C21") * sn(phi - g("phi21")))
                             * alpha ** 2
                             + 0.25 * (4 * g("C34") * sn(4 * (phi - g("phi34"))) + 2 * g("C32") * sn(2 * (phi - g("phi32"))))
                             * alpha ** 3
                             + 0.2 * (5 * g("C45") * sn(5 * (phi - g("phi45"))) + 3 * g("C43") * sn(3 * (phi - g("phi43")))
                                      + g("C41") * sn(phi - g("phi41"))) * alpha ** 4
                             + 1 / 6 * (6 * g("C56") * sn(6 * (phi - g("phi56"))) + 4 * g("C54") * sn(4 * (phi - g("phi54")))
                                        + 2 * g("C52") * sn(2 * (phi - g("phi52")))) * alpha ** 5)
    return np.exp(-np.sign(b) * (b / 2) ** 2 * (dk ** 2 + dp ** 2))


def text_temporal(alpha, lam, delta):  # 3.3
    return np.exp(-((0.5 * np.pi / lam * delta * alpha ** 2) ** 2))


def text_kernel(gpts, extent, energy, c, cutoff, soft, delta, beta):  # 3.5
    lam = s1.energy2wavelength(energy)
    sampling = (extent[0] / gpts[0], extent[1] / gpts[1])
    kx, ky = s1.fftfreq2(gpts, sampling)
    alpha = np.sqrt(kx ** 2 + ky ** 2) * lam
    phi = np.arctan2(ky + 0 * kx, kx + 0 * ky)
    h = np.exp(-1j * text_chi(alpha, phi, c, lam))
    if beta != 0:
        h = h * text_spatial(alpha, phi, c, lam, beta)
    if delta != 0:
        h = h * text_temporal(alpha, lam, delta)
    if np.isfinite(cutoff):
        if soft:
            dax = 1 / (gpts[0] * sampling[0]) * 1e3 * lam * 1e-3
            day = 1 / (gpts[1] * sampling[1]) * 1e3 * lam * 1e-3
            den = np.sqrt((np.cos(phi) * dax) ** 2 + (np.sin(phi) * day) ** 2)
            den[0, 0] = 1.0
            ap = np.clip((cutoff * 1e-3 - alpha) / den + 0.5, 0.0, 1.0)
            ap[0, 0] = 1.0
        else:
            ap = (alpha <= cutoff * 1e-3).astype(float)
        h = h * ap
    return h


def abtem_kernel(gpts, extent, energy, c, cutoff, soft, delta, beta):
    ctf = abtem.CTF(semiangle_cutoff=cutoff, soft=soft, focal_spread=delta, angular_spread=beta, energy=energy,
                    gpts=gpts, extent=extent, **c)
    return np.asarray(ctf._evaluate_kernel())


def main():
    pins = {"versions": {"abtem": abtem.__version__, "abtem_commit": "164e644ff411167db23ec4ef26aeb4fe8f384b26"},
            "cases": {}, "closed_forms": []}
    # Closed forms (3.7)
    for cs in (-5e5, 1e6, 1.2e7):
        for e in (80e3, 200e3, 300e3):
            lam = s1.energy2wavelength(e)
            sd, pr = float(scherzer_defocus(cs, e)), float(point_resolution(cs, e))
            check(f"scherzer {cs} {e}", abs(sd - np.sign(cs) * np.sqrt(1.5 * abs(cs) * lam)) / abs(sd), 1e-15)
            check(f"point resolution {cs} {e}", abs(pr - (lam ** 3 * abs(cs) / 6) ** 0.25) / pr, 1e-15)
            pins["closed_forms"].append({"Cs": cs, "energy": e, "scherzer": sd, "point_resolution": pr})
    pins["nyquist"] = []
    for a in (10.0, 20.0, 30.0):
        for e in (80e3, 200e3, 300e3):
            ns = float(nyquist_sampling(a, e))
            check(f"nyquist {a} {e}", abs(ns - 1 / (4 * a / s1.energy2wavelength(e) * 1e-3)) / ns, 1e-15)
            pins["nyquist"].append({"semiangle": a, "energy": e, "sampling": ns})

    # H1-H3: SrTiO3 exit wave, then the CTF and the image
    a0 = 3.905
    sto = Atoms("SrTiO3", scaled_positions=[(0, 0, 0), (0.5, 0.5, 0.5), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)],
                cell=[a0, a0, a0], pbc=True) * (3, 3, 5)
    gpts, energy = (120, 120), 300e3
    extent = (3 * a0, 3 * a0)
    pot = abtem.Potential(sto, gpts=gpts, slice_thickness=1.9525)
    exit_ = abtem.PlaneWave(energy=energy).multislice(pot).compute()
    psi = np.asarray(exit_.array)
    sch = float(scherzer_defocus(1.2e7, energy))
    cases = {
        "H1": ({"C30": 1.2e7, "C10": -sch}, 20.0, True, 0.0, 0.0),
        "H2": ({"C30": 1.2e7, "C10": -sch}, 20.0, True, 30.0, 0.5),
        "H3": ({"C10": -80.0, "C12": 50.0, "phi12": 0.4, "C21": 300.0, "phi21": 1.1, "C23": 200.0, "phi23": 0.1,
                "C30": 1.2e7}, np.inf, True, 40.0, 0.3),
    }
    for name, (c, cutoff, soft, delta, beta) in cases.items():
        ka = abtem_kernel(gpts, extent, energy, c, cutoff, soft, delta, beta)
        kt = text_kernel(gpts, extent, energy, c, cutoff, soft, delta, beta)
        check(f"{name} kernel text vs abTEM", s1.relmax(kt, ka), 1e-13)
        ctf = abtem.CTF(semiangle_cutoff=cutoff, soft=soft, focal_spread=delta, angular_spread=beta, **c)
        img = np.asarray(exit_.apply_ctf(ctf).intensity().compute().array)
        timg = np.abs(np.fft.ifft2(np.fft.fft2(psi) * kt)) ** 2
        check(f"{name} image text vs abTEM", s1.relmax(timg, img), 1e-12)
        put(f"{name}.kernel", ka)
        put(f"{name}.image", img)
        pins["cases"][name] = {"aberrations": c, "cutoff_mrad": None if not np.isfinite(cutoff) else cutoff, "soft": soft,
                               "focal_spread": delta, "angular_spread": beta, "gpts": list(gpts),
                               "extent": list(extent), "energy": energy, "image_from": "SrTiO3 3x3x5, slice 1.9525"}
    # K1, K2: kernels only
    allc = {"C10": 35.0, "C12": 12.0, "phi12": 0.3, "C21": 150.0, "phi21": -0.7, "C23": 90.0, "phi23": 1.9,
            "C30": -4e5, "C32": 2e4, "phi32": 0.6, "C34": 1e4, "phi34": -1.2, "C41": 3e6, "phi41": 0.2,
            "C43": 2e6, "phi43": 2.5, "C45": 1e6, "phi45": -0.4, "C50": 5e8, "C52": 2e8, "phi52": 0.9,
            "C54": 1e8, "phi54": -2.0, "C56": 5e7, "phi56": 0.05}
    kcases = {"K1": ((96, 160), (9.6, 14.4), 80e3, allc, 35.0, False, 25.0, 0.8),
              "K2": ((128, 128), (12.8, 12.8), 200e3, {"C10": -120.0, "C30": 1e6}, np.inf, True, 0.0, 0.4)}
    for name, (g, ext, e, c, cutoff, soft, delta, beta) in kcases.items():
        ka = abtem_kernel(g, ext, e, c, cutoff, soft, delta, beta)
        kt = text_kernel(g, ext, e, c, cutoff, soft, delta, beta)
        check(f"{name} kernel text vs abTEM", s1.relmax(kt, ka), 1e-13)
        put(f"{name}.kernel", ka)
        pins["cases"][name] = {"aberrations": c, "cutoff_mrad": None if not np.isfinite(cutoff) else cutoff, "soft": soft,
                               "focal_spread": delta, "angular_spread": beta, "gpts": list(g), "extent": list(ext),
                               "energy": e}

    os.makedirs(FIXTURES, exist_ok=True)
    blob = bytes(BLOB)
    open(os.path.join(FIXTURES, "slice-3-arrays.bin"), "wb").write(blob)
    pins["arrays"] = {"file": "slice-3-arrays.bin", "sha256": hashlib.sha256(blob).hexdigest(), "bytes": len(blob),
                      "index": INDEX}
    pins["self_checks"] = CHECKS
    json.dump(pins, open(os.path.join(FIXTURES, "slice-3-pins.json"), "w"), indent=1, sort_keys=True)
    print("\n".join(CHECKS))
    print(f"wrote {len(INDEX)} arrays, {len(blob) / 1e6:.2f} MB; {len(CHECKS)} self-checks passed")


if __name__ == "__main__":
    main()
