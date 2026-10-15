"""SimulatoKit performance lane, part 2: the break-once proof (registration-performance.md section 6, gate 3).

The runner is tools/slice-1/mutants.py; this file holds the lane's mutants. Usage:  python3 tools/performance/mutants.py <clone> [ids...]
"""

import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m1", os.path.join(HERE, "..", "slice-1", "mutants.py"))
m1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m1)

E = "Sources/SimulatoMultislice/MetalEngine.swift"
POT = "Sources/SimulatoMultislice/Potential.swift"
PROBE = "Sources/SimulatoProbe/Probe.swift"
F = "MetalGPUTests|MetalTests|STEMPinsTests|PotentialTests"

m1.MUTANTS = [
    ("P01", E, "let phase = g.multiplication(v, sigma, name: nil)",
     "let phase = g.multiplication(v, g.constant(1e-3, dataType: .float32), name: nil)", F, "σ replaced in the GPU phase"),
    ("P02", E, "let kernel = try builder.kernel(z)", "let kernel = try builder.kernel(elements[(q + 1) % elements.count])",
     F, "the kernel of the next element"),
    ("P03", E, "let w = [1 + xy - y - x, x - xy, y - xy, xy]", "let w = [1 + xy - y - x, y - xy, x - xy, xy]", F,
     "the delta weights' x and y swapped"),
    ("P04", E, "if let p = byThickness[dz] { return p }", "if let p = byThickness.values.first { return p }", F,
     "the first thickness's propagator for every slice"),
    ("P05", PROBE, "let sx = -2 * Double.pi * px[ix] * ux\n            rx.real[ix]",
     "let sx = 2 * Double.pi * px[ix] * ux\n            rx.real[ix]", F, "the GPU probe's x ramp sign flipped"),
    ("P06", E, "trueTensor: g.squareRoot(with: g.division(num, den, name: nil), name: nil),",
     "trueTensor: g.division(num, den, name: nil),", F, "the renormalisation squared (s² for s)"),
    ("P07", E, "let i2 = g.addition(g.multiplication(re, re, name: nil), g.multiplication(im, im, name: nil), name: nil)",
     "let i2 = g.multiplication(re, re, name: nil)", F, "the scan intensity as Re² only"),
    ("P08", E, "                try read(j, current)\n", "                try read(j, other)\n", F,
     "the wave read one slice early"),
    ("P09", E, "perSlice[j][e * n + i * ny + jj, default: 0] += w[k]", "perSlice[j][e * n + i * ny + jj] = w[k]", F,
     "the host merge overwrites instead of summing"),
    ("P10", PROBE, "            }\n        }\n        let f = norm.squareRoot()\n        for i in psi.real.indices {\n            psi.real[i] /= f\n            psi.imaginary[i] /= f\n        }\n        return psi",
     "            }\n        }\n        let f = norm.squareRoot() * 1.0001\n        for i in psi.real.indices {\n            psi.real[i] /= f\n            psi.imaginary[i] /= f\n        }\n        return psi",
     F, "the GPU probe's spectrum normalised 1e−4 off"),
    ("P11", E, "let sum = g.reductionSum(with: fk, axis: 0, name: nil)",
     "let sum = g.sliceTensor(fk, dimension: 0, start: 0, length: 1, name: nil)", F, "only the first element summed"),
    ("P12", E, "let limited = g.multiplication(spectrum, g.cast(aperture, to: .complexFloat32, name: \"a\"), name: nil)",
     "let limited = g.multiplication(spectrum, g.constant(1.0, dataType: .complexFloat32), name: nil)", F,
     "the GPU transmission not bandlimited"),
    ("P13", POT, "values: try builder.potential(of: atoms(inSlice: j))", "values: try builder.potential(of: atoms(inSlice: 0))",
     F, "Potential.build puts slice 0's atoms in every slice"),
    ("P14", E, "var idx = [Int32](repeating: 0, count: m), w = [Float](repeating: 0, count: m)",
     "var idx = [Int32](repeating: 0, count: m), w = [Float](repeating: 1e-3, count: m)", F,
     "the padding carries weight"),
    ("P15", "Sources/SimulatoMultislice/ProjectedPotential.swift",
     "                             parameterization: parameterization).build()",
     "                             parameterization: LobatoParameterization()).build()", "PotentialTests",
     "the eager potential drops the parameterisation it was given"),
    # Gate B's fixes (amendment 3).
    ("P16", E, "        guard p.sliceCount > 0 else { return try read(-1, psi) }\n", "", "MetalGPUTests",
     "the no-slice guard removed"),
    ("P17", E, "let ratio = g.multiplication(g.division(eY, eX, name: nil), g.division(ex, ey, name: nil), name: nil)",
     "let ratio = g.division(g.multiplication(eY, ex, name: nil), g.multiplication(eX, ey, name: nil), name: nil)",
     "MetalGPUTests", "the renormalisation by products again"),
    ("P18", E, "    static let cacheLimit = 12", "    static let cacheLimit = 100_000", "MetalGPUTests",
     "the executable cache unbounded"),
    # Amendment 4: the accepted floor of the detector bound.
    ("P19", "Sources/SimulatoMultislice/ExecutionOptions.swift", "                detectorFloor: Double = 1e-5,",
     "                detectorFloor: Double = 0,", "STEMPinsTests|GeneralityTests", "the detector floor at 0 (the old bound)"),
]
m1.EXTRA = {}

if __name__ == "__main__":
    m1.main()
