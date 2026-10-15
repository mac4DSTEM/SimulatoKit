"""SimulatoKit slice 1: the break-once proof (registration section 6, gate 4).

Each mutant is one textual change aimed at one registered rule, applied to a clone of the frozen commit; the tests named
for it must fail (a compile error does not count). Usage, from any directory, with a clone of the commit to mutate:

    python3 tools/slice-1/mutants.py <clone> > mutants.log

The clone is reset (`git checkout -- .`) before every mutant. One mutant at a time, `swift test -c release --filter`.
Prints one line per mutant: KILLED (the named tests failed), SURVIVED, or BUILD-FAILED, and the tests that failed (Swift
Testing's "✘ Test name(" lines).
"""

import re
import subprocess
import sys

MUTANTS = [
    # (id, file, old, new, filter, intent)
    ("M01", "Sources/SimulatoCore/ElectronBeam.swift",
     "(energy * (2 * c.electronMass * c2 / c.elementaryCharge + energy))",
     "(energy * (2 * c.electronMass * c2 / c.elementaryCharge))", "BeamAndScatteringTests",
     "λ without the relativistic term"),
    ("M02", "Sources/SimulatoCore/ElectronBeam.swift",
     "2 * Double.pi * mass * ASEUnits.kg", "2 * Double.pi * c.electronMass * ASEUnits.kg", "BeamAndScatteringTests",
     "σ with the rest mass"),
    ("M03", "Sources/SimulatoMultislice/LobatoParameterization.swift",
     "        return (a, b)\n", "        return (b, a)\n", "BeamAndScatteringTests|PotentialPinsTests",
     "A and B swapped"),
    ("M04", "Sources/SimulatoMultislice/ProjectedPotential.swift",
     "out[i] = f[i] / compensation[i]", "out[i] = f[i] / (grid.samplingX * grid.samplingY)", "PotentialPinsTests|TruthTests",
     "the compensation S dropped"),
    ("M05", "Sources/SimulatoMultislice/ProjectedPotential.swift",
     "let k = (a * a + b * b).squareRoot()", "let k = Double.pi * (a * a + b * b).squareRoot()",
     "PotentialPinsTests|TruthTests", "S with π inside"),
    ("M06", "Sources/SimulatoMultislice/ProjectedPotential.swift",
     "w[1][a] = x - xy\n            w[2][a] = y - xy", "w[1][a] = y - xy\n            w[2][a] = x - xy",
     "PotentialPinsTests", "the bilinear weights' fx and fy swapped"),
    ("M07", "Sources/SimulatoMultislice/Multislice.swift",
     "let x = -k2 * Double.pi * dz * lambda", "let x = k2 * Double.pi * dz * lambda",
     "MultislicePinsTests|TruthTests", "the propagator's sign flipped"),
    ("M08", "Sources/SimulatoMultislice/Multislice.swift",
     "propagators[slice.thickness] = propagator(grid: grid, beam: incident.beam, thickness: slice.thickness)",
     "propagators[slice.thickness] = propagator(grid: grid, beam: incident.beam, thickness: potential.thickness)",
     "MultislicePinsTests|TruthTests", "the propagator over the whole thickness"),
    ("M09", "Sources/SimulatoCore/Antialias.swift",
     "public static let cutoff = 0.6666666", "public static let cutoff = 0.5",
     "PotentialPinsTests|MultislicePinsTests", "the antialias cutoff 1/2"),
    ("M10", "Sources/SimulatoMultislice/Multislice.swift",
     "                t.real[i] *= aperture[i]\n                t.imaginary[i] *= aperture[i]\n",
     "                _ = aperture[i]\n", "PotentialPinsTests|MultislicePinsTests", "the transmission not bandlimited"),
    ("M11", "Sources/SimulatoMultislice/Multislice.swift",
     "            multiply(&psi, by: t)\n            fft.forward(&psi)",
     "            fft.forward(&psi)", "MultislicePinsTests",
     "propagate before transmit (part 1: transmit removed here, re-added after the inverse by M11b's second edit)"),
    ("M12", "Sources/SimulatoMultislice/Slicing.swift",
     "if p.z > box.z - 1e-10 { p.z = 0 }", "if p.z > box.z + 1 { p.z = 0 }", "PotentialPinsTests",
     "the top-face rule dropped"),
    ("M13", "Sources/SimulatoMultislice/Slicing.swift",
     "edges = edges.map { $0 - 1e-12 }", "edges = edges.map { $0 }", "PotentialPinsTests",
     "the slice-edge nudge dropped"),
    ("M14", "Sources/SimulatoProbe/Probe.swift",
     "public init(semiangleCutoff: Double, soft: Bool = true)", "public init(semiangleCutoff: Double, soft: Bool = false)",
     "MultislicePinsTests", "the probe aperture hard"),
    ("M15", "Sources/SimulatoProbe/Probe.swift",
     "        let f = norm.squareRoot()\n        for i in psi.real.indices {\n            psi.real[i] /= f\n            psi.imaginary[i] /= f\n        }\n        try FFT2D",
     "        let f = norm.squareRoot() / Double(grid.count).squareRoot()\n        for i in psi.real.indices {\n            psi.real[i] /= f\n            psi.imaginary[i] /= f\n        }\n        try FFT2D",
     "MultislicePinsTests|TruthTests", "the probe normalised in real space"),
    ("M16", "Sources/SimulatoProbe/Probe.swift",
     "let (c, s) = (cos(-chi), sin(-chi))", "let (c, s) = (cos(chi), sin(chi))", "MultislicePinsTests",
     "χ's sign flipped"),
    ("M17", "Sources/SimulatoProbe/Probe.swift",
     "let sx = -2 * Double.pi * px[ix] * ux\n            let (cx", "let sx = 2 * Double.pi * px[ix] * ux\n            let (cx", "MultislicePinsTests|TruthTests",
     "the shift kernel's sign flipped"),
    ("M18", "Sources/SimulatoCore/Antialias.swift",
     "            if m % 2 != n % 2 { m += 1 }\n", "", "MultislicePinsTests|SimulatoCoreTests",
     "the crop's first parity step dropped"),
    ("M19", "Sources/SimulatoCore/Waves.swift",
     "if normalization == .values {", "if normalization == .reciprocalSpace {", "MultislicePinsTests|TruthTests",
     "the plane wave's 1/N dropped in the pattern"),
    ("M20", "Sources/SimulatoMultislice/MetalEngine.swift",
     "let propagated = g.multiplication(spectrum, p, name: nil)",
     "let propagated = g.multiplication(spectrum, g.constant(1.0, dataType: .complexFloat32), name: nil)",
     "MetalTests", "Metal: the propagator not applied"),
    ("M21", "Sources/SimulatoMultislice/MetalEngine.swift",
     "        d.scalingMode = .size", "        d.scalingMode = .unitary", "MetalTests",
     "Metal: the inverse FFT scaled by 1/sqrt(N)"),
    ("M22", "Sources/SimulatoMultislice/LobatoParameterization.swift",
     "case .abTEMFloat32: (pi, pi2) = (Self.float32Pi, Self.float32PiSquared)",
     "case .abTEMFloat32: (pi, pi2) = (Double.pi, Double.pi * Double.pi)", "PotentialPinsTests|MultislicePinsTests",
     "abTEM's float32 π replaced by the exact π"),
    ("M23", "Sources/SimulatoMultislice/LobatoParameterization.swift",
     "public static let float32Pi = Double(Float(Double.pi))", "public static let float32Pi = Double(Float.pi)",
     "BeamAndScatteringTests", "Swift's Float.pi (rounded toward zero) for the float32 π"),
    ("M24", "Sources/SimulatoMultislice/LobatoParameterization.swift",
     "case .exact: (pi, pi2) = (Double.pi, Double.pi * Double.pi)",
     "case .exact: (pi, pi2) = (Self.float32Pi, Self.float32PiSquared)", "PotentialPinsTests|BeamAndScatteringTests",
     "the exact π replaced by abTEM's"),
    ("M25", "Sources/SimulatoCore/ElectronBeam.swift",
     "guard energy.isFinite, energy > 0 else", "guard energy.isFinite, energy >= 0 else", "BeamAndScatteringTests",
     "a zero energy accepted"),
    ("M26", "Sources/SimulatoMultislice/ProjectedPotential.swift",
     "            w[3][a] = xy\n", "            w[3][a] = 0.5 * xy\n", "PotentialPinsTests",
     "the bilinear weights no longer sum to one"),
    ("M27", "Sources/SimulatoMultislice/Slicing.swift",
     "let n = (depth / s).rounded(.up)", "let n = (depth / s).rounded(.down)", "PotentialPinsTests",
     "⌊Lz / s⌋ slices instead of ⌈Lz / s⌉"),
    ("M28", "Sources/SimulatoMultislice/Slicing.swift",
     "guard abs(sum - depth) <= 1e-8 + 1e-5 * abs(depth) else", "guard abs(sum - depth) <= 1e3 else",
     "PotentialPinsTests", "slice thicknesses that do not sum to the depth accepted"),
    ("M29", "Sources/SimulatoCore/Antialias.swift",
     "                    value = 1\n", "                    value = 0.999999\n", "TruthTests",
     "the aperture below one inside the cutoff"),
    ("M30", "Sources/SimulatoMultislice/Multislice.swift",
     "let k2 = kx[ix] * kx[ix] + ky[iy] * ky[iy]", "let k2 = kx[ix] * kx[ix] + 2 * ky[iy] * ky[iy]", "TruthTests",
     "an anisotropic propagator (breaks 4mm)"),
    ("M31", "Sources/SimulatoMultislice/ProjectedPotential.swift",
     "SIMD2($0.position.x / grid.samplingX,", "SIMD2(($0.position.x + 0.1 * $0.position.z) / grid.samplingX,",
     "TruthTests", "atoms shifted along x by depth (breaks the slice lattice)"),
    ("M32", "Sources/SimulatoMultislice/Multislice.swift",
     "p.real[i] = cos(x) * aperture[i]", "p.real[i] = cos(x) * aperture[i] * 1.001", "TruthTests",
     "a propagator with gain (Σ|ψ|² grows)"),
    ("M33", "Sources/SimulatoProbe/Probe.swift",
     "let sy = -2 * Double.pi * py[iy] * uy\n                let (cy", "let sy = -2 * Double.pi * py[iy] * (uy + 0.37)\n                let (cy", "TruthTests",
     "the probe off its position along y"),
    ("M34", "Sources/SimulatoProbe/Aberrations.swift",
     "x = x + 1 / 3 * pow(alpha, 3)", "x = x + 1 / 3 * pow(alpha, 2)", "TruthTests",
     "χ's second-order term with α²"),
    ("M35", "Sources/SimulatoCore/Grid.swift",
     "return [1, 3, 5, 15].contains(m) && power >= 3", "return [1, 3, 5, 15].contains(m) && power >= 2",
     "SimulatoCoreTests", "lengths vDSP refuses accepted"),
    ("M36", "Sources/SimulatoCore/Grid.swift",
     "let positive = (n - 1) / 2 + 1", "let positive = n / 2 + 1", "SimulatoCoreTests",
     "fftfreq's Nyquist bin positive"),
    ("M37", "Sources/SimulatoCore/FFT2D.swift",
     "transform(&field.real, &field.imaginary, rows: forwardY, columns: forwardX)",
     "transform(&field.real, &field.imaginary, rows: inverseY, columns: forwardX)", "SimulatoCoreTests",
     "the forward FFT inverse along y"),
    ("M38", "Sources/SimulatoCore/AtomicStructure.swift",
     "let shift = SIMD3<Double>(Double(mx) * box.x, Double(my) * box.y, Double(mz) * box.z)",
     "let shift = SIMD3<Double>(Double(mz) * box.x, Double(my) * box.y, Double(mx) * box.z)", "SimulatoCoreTests",
     "the repeat's axes swapped"),
    ("M39", "Sources/SimulatoProbe/Probe.swift",
     "guard semiangleCutoff > 0, !semiangleCutoff.isNaN else", "guard !semiangleCutoff.isNaN else", "SimulatoProbeTests",
     "a zero semiangle accepted"),
    ("M40", "Sources/SimulatoDetector/SimulatoDetector.swift",
     "        \"bright-field, annular (ADF, HAADF), pixelated (4D) and segmented detectors with frozen phonons (slice 2)\"",
     "        \"\"", "SimulatoDetectorPlaceholderTests", "the placeholder's scope emptied"),
]

# M11 needs a second edit (transmit after the inverse FFT); applied together.
EXTRA = {
    "M11": ("Sources/SimulatoMultislice/Multislice.swift",
            "            fft.inverse(&psi)\n            if wanted.contains(j) {",
            "            fft.inverse(&psi)\n            multiply(&psi, by: t)\n            if wanted.contains(j) {"),
}


def sh(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, shell=True, capture_output=True, text=True)


def apply(clone, path, old, new):
    full = f"{clone}/{path}"
    text = open(full).read()
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"{path}: the mutant's text occurs {n} times, not once: {old!r}")
    open(full, "w").write(text.replace(old, new))


def main():
    clone = sys.argv[1]
    only = set(sys.argv[2:])
    for mid, path, old, new, flt, intent in MUTANTS:
        if only and mid not in only:
            continue
        sh("git checkout -- .", clone)
        apply(clone, path, old, new)
        if mid in EXTRA:
            apply(clone, *EXTRA[mid])
        r = sh(f"swift test -c release --filter '{flt}' 2>&1", clone)
        out = r.stdout
        # Swift Testing (the kits' standard since 2026-10-15): "✘ Test name(args) ... failed" / "recorded an issue".
        failed = sorted(set(re.findall(r"✘ Test (\w+)\(", out)))
        if "Test run started" in out and r.returncode != 0 and not failed:
            verdict = "KILLED"  # the test process crashed after the run started: the tests did not pass
            failed = ["(test process crashed)"]
        elif "error: compile" in out or ("error:" in out and "Test run with" not in out and not failed):
            verdict = "BUILD-FAILED"
        elif r.returncode != 0 and failed:
            verdict = "KILLED"
        elif r.returncode != 0:
            verdict = "BUILD-FAILED"
        else:
            verdict = "SURVIVED"
        print(f"{mid} {verdict:12} {intent} [{flt}] failed: {', '.join(failed) if failed else '-'}", flush=True)
    sh("git checkout -- .", clone)


if __name__ == "__main__":
    main()
