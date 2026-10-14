"""SimulatoKit slice 2: the break-once proof (registration-slice-2.md section 6, gate 4).

The runner is tools/slice-1/mutants.py (one textual change at a time on a clone, `swift test -c release --filter`, Swift
Testing's output); this file holds slice 2's mutants. Usage:  python3 tools/slice-2/mutants.py <clone> [ids...]
"""

import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m1", os.path.join(HERE, "..", "slice-1", "mutants.py"))
m1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m1)

FP = "Sources/SimulatoMultislice/FrozenPhonons.swift"
DET = "Sources/SimulatoDetector/Detectors.swift"
STEM = "Sources/SimulatoDetector/STEM.swift"
SCAN = "Sources/SimulatoProbe/GridScan.swift"
MS = "Sources/SimulatoMultislice/MetalMultislice.swift"

m1.MUTANTS = [
    ("N01", FP, "position: $0.0.position + $0.1)", "position: $0.0.position + 2 * $0.1)", "FrozenPhononsTests|STEMPinsTests",
     "σ applied twice (displacement doubled)"),
    ("N02", FP, "return SIMD3(directions.contains(.x) ? s * r.x : 0,", "return SIMD3(s * r.x,", "FrozenPhononsTests",
     "the directions ignored on x"),
    ("N03", FP, "var g = GaussianGenerator(seed: seed &+ UInt64(c))", "var g = GaussianGenerator(seed: seed)",
     "FrozenPhononsTests", "every configuration from the same seed"),
    ("N04", FP, "let u2 = Double(nextUInt64() >> 11) * 0x1p-53", "let u2 = Double(nextUInt64() >> 11) * 0x1p-54",
     "FrozenPhononsTests", "Box–Muller's angle over half the circle (wrong statistics)"),
    ("N05", FP, "spare = r * sin(theta)", "spare = r * cos(theta)", "FrozenPhononsTests",
     "the spare normal equal to the first (correlated pairs)"),
    ("N06", SCAN, "let div = endpoint ? n - 1 : n", "let div = endpoint ? n : n - 1", "STEMPinsTests",
     "linspace with endpoint inverted"),
    ("N07", DET, "if a2 >= i2 && a2 < o2 { sum += p.values[ix * p.grid.ny + iy] }",
     "if a2 >= i2 && a2 <= o2 { sum += p.values[ix * p.grid.ny + iy] }", "DetectorTruthTests|STEMPinsTests",
     "the annular mask closed at outer"),
    ("N08", DET, "if a2 >= i2 && a2 < o2 { sum += p.values[ix * p.grid.ny + iy] }",
     "if a2 > i2 && a2 < o2 { sum += p.values[ix * p.grid.ny + iy] }", "DetectorTruthTests|STEMPinsTests",
     "inner as > instead of ≥"),
    ("N09", DET, "        if waves.normalization == .values {\n            let n = Double(waves.grid.count)",
     "        if waves.normalization == .reciprocalSpace {\n            let n = Double(waves.grid.count)",
     "STEMPinsTests|DetectorTruthTests", "the probe's pattern divided by N"),
    ("N10", DET, "guard let o = outer else { let c = p.cutoffAngles; return floor(min(c.x, c.y)) }",
     "guard let o = outer else { let c = p.cutoffAngles; return min(c.x, c.y) }", "STEMPinsTests|DetectorTruthTests",
     "outer = nil without the floor"),
    ("N11", DET, "var k = Int(floor(DiffractionIntensity.roundToDecimals(Double(m) / 2.0.squareRoot(), 7)))",
     "var k = Int(floor(DiffractionIntensity.roundToDecimals(Double(m) / 1.5, 7)))", "STEMPinsTests",
     "pixelated 'valid' without √2"),
    ("N12", DET, "var phi = (atan2(ay[cy], ax[cx]) - rotation)", "var phi = (atan2(ax[cx], ay[cy]) - rotation)",
     "STEMPinsTests", "the segment's φ from atan2(kx, ky)"),
    ("N13", DET, "var phi = (atan2(ay[cy], ax[cx]) - rotation)", "var phi = (atan2(ay[cy], ax[cx]) + rotation)",
     "STEMPinsTests", "the rotation added instead of subtracted"),
    ("N14", DET, "let r = Int((Double(radialBins) * (alpha - inner) / (outer - inner)).rounded(.down))",
     "let r = min(radialBins - 1, Int((Double(radialBins) * (alpha - inner) / (outer - inner)).rounded()))",
     "STEMPinsTests", "the radial bin rounded instead of floored"),
    ("N15", STEM, "values: values.map { $0 / n })", "values: values)", "STEMPinsTests",
     "the average over configurations a sum"),
    ("N16", MS, "let spectrum = graph.fastFourierTransform(transmitted, axes: [1, 2], descriptor: forward, name: nil)",
     "let spectrum = graph.fastFourierTransform(transmitted, axes: [0, 1], descriptor: forward, name: nil)",
     "STEMPinsTests", "Metal: the batch's forward FFT over axes 0 and 1"),
    ("N17", DET, "        if o > maximum && abs(o - maximum) > 1e-5 + 1e-8 * abs(maximum) {",
     "        if o > 10 * maximum {", "DetectorTruthTests", "an outer beyond the maximum angle accepted"),
    ("N18", DET, "            if m % 2 != n % 2 { m += 1 }\n            return m\n",
     "            return m\n", "STEMPinsTests", "the cutoff gpts without the grid's parity"),
    ("N19", DET, "        let mx = m(p.grid.nx, a.x), my = m(p.grid.ny, a.y)",
     "        let mx = p.grid.nx, my = p.grid.ny", "STEMPinsTests", "the segmented pattern not cropped at outer"),
    ("N20", "Sources/SimulatoDetector/SimulatoDetector.swift",
     "        \"bright-field, annular (ADF, HAADF), pixelated (4D) and segmented detectors with frozen phonons (slice 2)\"",
     "        \"\"", "DetectorTruthTests", "the module's scope emptied"),
    ("N21", "Sources/SimulatoMultislice/ProjectedPotential.swift", "            let k = try kernel(z)\n",
     "            let k = try kernel(6)\n", "DetectorTruthTests", "every element scattered as carbon"),
    ("N22", STEM, "for k in 0..<v.count { sums[d][q * counts[d] + k] += v[k] }",
     "for k in 0..<v.count { sums[d][q * counts[d] + k] += v[k] * Double(q + 1) }", "DetectorTruthTests",
     "the scan weighted by position"),
    ("N23", FP, "position: $0.0.position + $0.1)", "position: $0.0.position + $0.1 + SIMD3(1e-9, 0, 0))",
     "FrozenPhononsTests", "every configuration offset by 1e-9 Å"),
]
m1.EXTRA = {}

if __name__ == "__main__":
    m1.main()
