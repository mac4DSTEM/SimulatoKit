"""SimulatoKit slice 3: the break-once proof (registration-slice-3.md section 6, gate 3).

The runner is tools/slice-1/mutants.py; this file holds slice 3's mutants. Usage:  python3 tools/slice-3/mutants.py <clone> [ids...]
"""

import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m1", os.path.join(HERE, "..", "slice-1", "mutants.py"))
m1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m1)

CTF = "Sources/SimulatoProbe/ContrastTransferFunction.swift"
F = "CTFTests|CTFTruthTests"

m1.MUTANTS = [
    ("C01", CTF, "        return exp(-(x * x))\n", "        return exp(-x)\n", F, "E_t without the square"),
    ("C02", CTF, "let x = 0.5 * Double.pi / wavelength * focalSpread * alpha * alpha",
     "let x = 0.5 * Double.pi * wavelength * focalSpread * alpha * alpha", F, "E_t with λ multiplied"),
    ("C03", CTF, "let b = angularSpread / 1e3", "let b = angularSpread", F, "β not converted from mrad"),
    ("C04", CTF, "let t1 = 1 / 2.0 * (2.0 * C12 * sin(2.0 * (p - phi12))) * a",
     "let t1 = (2.0 * C12 * sin(2.0 * (p - phi12))) * a", F, "∂χ/∂φ's first ½ dropped"),
    ("C05", CTF, "let t3 = (C34 * cos(4 * (p - phi34)) + C32 * cos(2 * (p - phi32)) + C30) * pow(a, 3)",
     "let t3 = (C34 * cos(4 * (p - phi34)) + C32 * cos(2 * (p - phi32))) * pow(a, 3)", F, "∂χ/∂α missing C30"),
    ("C06", CTF, "        return exp(-sign * (b / 2) * (b / 2) * (dk * dk + dp * dp))",
     "        return exp(-(b / 2) * (b / 2) * (dk * dk + dp * dp)) + 0 * sign", F, "sign(β) dropped (equivalent for β ≥ 0)"),
    ("C07", CTF, "? try Aperture(semiangleCutoff: semiangleCutoff, soft: softAperture)",
     "? try Aperture(semiangleCutoff: semiangleCutoff, soft: false)", F, "the soft aperture made hard"),
    ("C08", CTF, "                var re = cos(-chi), im = sin(-chi)", "                var re = cos(chi), im = sin(chi)", F,
     "the kernel conjugated"),
    ("C09", CTF, "                    factor *= Self.temporalEnvelope(alpha: alpha, wavelength: wavelength, focalSpread: focalSpread)\n",
     "                    factor *= Self.temporalEnvelope(alpha: alpha, wavelength: wavelength, focalSpread: focalSpread)\n"
     "                    factor *= Self.temporalEnvelope(alpha: alpha, wavelength: wavelength, focalSpread: focalSpread)\n",
     F, "the temporal envelope applied twice"),
    ("C10", CTF, "        fft.inverse(&f)\n        return Waves(field: f", "        fft.inverse(&f)\n        f.real = f.real.map { $0 * 2 }\n        return Waves(field: f",
     F, "apply_ctf scaled wrongly (×2 on the real part)"),
    ("C11", CTF, "return s * (3 / 2 * abs(Cs) * lambda).squareRoot()", "return s * (4 / 3 * abs(Cs) * lambda).squareRoot()", F,
     "Scherzer's 3/2 replaced by 4/3"),
    ("C12", CTF, "return pow(pow(lambda, 3) * abs(Cs) / 6, 1.0 / 4)", "return pow(pow(lambda, 3) * abs(Cs) / 6, 1.0 / 3)", F,
     "point resolution with a cube root"),
    ("C13", CTF, "return 1 / (4 * semiangleCutoff / lambda * 1e-3)", "return 1 / (2 * semiangleCutoff / lambda * 1e-3)", F,
     "Nyquist with 2 instead of 4"),
    ("C14", CTF, "guard angularSpread.isFinite, angularSpread >= 0 else", "guard angularSpread.isFinite else", F,
     "a negative angular spread accepted"),
    ("C15", CTF, "                if let a = aperture { factor *= a[i] }\n", "", F, "the aperture not applied"),
]
m1.EXTRA = {}

if __name__ == "__main__":
    m1.main()
