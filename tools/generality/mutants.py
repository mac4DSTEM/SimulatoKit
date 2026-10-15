"""SimulatoKit 0.2.0, the generality pass: the break-once proof (registration-generality.md section 3).

The runner is tools/slice-1/mutants.py; this file holds the pass's mutants. Usage:  python3 tools/generality/mutants.py <clone> [ids...]
"""

import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("m1", os.path.join(HERE, "..", "slice-1", "mutants.py"))
m1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m1)

E = "Sources/SimulatoMultislice/MetalEngine.swift"
O = "Sources/SimulatoMultislice/ExecutionOptions.swift"
CTF = "Sources/SimulatoProbe/ContrastTransferFunction.swift"
F = "GeneralityTests|CTFConventionTests"

m1.MUTANTS = [
    ("G01", E, "let size = batchSize(grid, bytes: options.metalBatchBytes)", "let size = batchSize(grid)", F,
     "the Metal batch bytes ignored"),
    ("G02", E, "        guard shapes == .powersOfTwo else { return b }\n", "", F, "the exact batch shapes ignored"),
    ("G03", E, "                           limit: options.metalExecutableCacheLimit) {", "                           limit: 12) {", F,
     "the executable cache limit ignored"),
    ("G04", E, "let out = MetalGraphs.propagate(g, psi: psi, t: t, p: prop, renormalize: options.metalRenormalization)",
     "let out = MetalGraphs.propagate(g, psi: psi, t: t, p: prop, renormalize: true)", F, "the renormalisation forced on"),
    ("G05", "Sources/SimulatoCore/Grid.swift", "abs(extentX - other.extentX) <= relativeTolerance *",
     "abs(extentX - other.extentX) <= 1e-12 *", F, "the grid tolerance ignored"),
    ("G06", "Sources/SimulatoMultislice/Potential.swift",
     "let thicknesses = try Slicing.validated(sliceThicknesses, depth: structure.box.z, maximumCount: maximumSliceCount)",
     "let thicknesses = try Slicing.validated(sliceThicknesses, depth: structure.box.z)", F, "the slice ceiling ignored"),
    ("G07", CTF, "return s * (factor * abs(Cs) * lambda).squareRoot()", "return s * (1.5 * abs(Cs) * lambda).squareRoot()", F,
     "the Scherzer factor ignored"),
    ("G08", CTF, "return pow(pow(lambda, 3) * abs(Cs) / divisor, 1.0 / 4)", "return pow(pow(lambda, 3) * abs(Cs) / 6, 1.0 / 4)", F,
     "the point-resolution divisor ignored"),
    ("G09", CTF, "return 1 / (oversampling * semiangleCutoff / lambda * 1e-3)", "return 1 / (4 * semiangleCutoff / lambda * 1e-3)",
     F, "the Nyquist oversampling ignored"),
    ("G10", O, "    public init(wave: Double = 1e-4,", "    public init(wave: Double = 2e-4,", F,
     "a registered Metal bound moved"),
    ("G11", O, "        self.maximumCPUThreads = max(1, maximumCPUThreads)", "        self.maximumCPUThreads = maximumCPUThreads", F,
     "the thread count no longer clamped"),
    # The tests that an option changes no number (or rounding only), broken by making the option change one.
    ("G12", "Sources/SimulatoMultislice/MetalMultislice.swift", "                    results.withLock { $0[q] = w }",
     "                    results.withLock { $0[q] = workers == 1 && q == count - 1 ? incidents[q] : w }", F,
     "one thread returns the last incident unpropagated"),
    ("G13", "Sources/SimulatoDetector/STEM.swift", "            let range = first..<min(first + chunk, positions.count)\n",
     "            let range = first..<min(first + chunk, positions.count)\n            if chunk == 1 && first == 0 { continue }\n",
     "theScanChunkIsAParameter", "a one-position chunk skips the first position"),
    ("G14", E, "out += fields(data, batch: b, grid: p.grid, count: chunk.count)",
     "out += fields(data, batch: b, grid: p.grid, count: b == used ? used : used - 1)", F,
     "a padded batch read one wave short"),
]
# Equivalent by design, not run: the thread cap and the CPU scan chunk change no number (bit for bit; tested as such).
m1.EXTRA = {}

if __name__ == "__main__":
    m1.main()
