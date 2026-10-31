//
//  EnsembleTests.swift
//  Role: the two ensemble axes abTEM adds to a scan: Poisson noise with `samples > 1` (a sample axis) and the exit waves
//        of every frozen-phonon configuration (a configuration axis, never averaged).
//

import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoKitTestSupport
import SimulatoMultislice
import SimulatoProbe

@Suite("Ensembles")
struct EnsembleTests {
    /// A 2 x 3 scan of one value per position.
    static func measurement() -> ScanMeasurement {
        ScanMeasurement(scanShape: [2, 3], valueShape: [], values: [0.5, 1, 2, 4, 8, 0])
    }

    @Test func poissonSamplesAreIndependentDrawsOfTheSameMean() throws {
        let flat = ScanMeasurement(scanShape: [20, 25], valueShape: [], values: [Double](repeating: 1, count: 500))
        let samples = try flat.poissonNoise(dose: 4, samples: 40, seed: 3)
        #expect(samples.count == 40)
        for s in samples { #expect(s.scanShape == [20, 25] && s.valueShape == [] && s.values.count == 500) }
        // 20 000 draws of mean 4: mean and variance 4 within about 5 standard errors (0.0141 and 0.0424).
        let all: [Double] = samples.flatMap { $0.values }
        let mean = all.reduce(0, +) / Double(all.count)
        let squares: [Double] = all.map { ($0 - mean) * ($0 - mean) }
        let variance = squares.reduce(0, +) / Double(all.count - 1)
        #expect(abs(mean - 4) < 0.07)
        #expect(abs(variance - 4) < 0.2)
        // The samples differ from one another, and the first is the single draw of the same seed.
        #expect(samples[0].values != samples[1].values)
        #expect(samples[0].values == (try flat.poissonNoise(dose: 4, seed: 3)).values)
        // Per sample at one position: the same mean across samples.
        let atZero: [Double] = samples.map { $0.values[0] }
        let meanAtZero = atZero.reduce(0, +) / 40
        #expect(abs(meanAtZero - 4) < 1.6)
    }

    @Test func poissonSamplesAreSeededAndRefuseBadCounts() throws {
        let m = Self.measurement()
        let a = try m.poissonNoise(dose: 3, samples: 3, seed: 9), b = try m.poissonNoise(dose: 3, samples: 3, seed: 9)
        #expect(a.map(\.values) == b.map(\.values))
        // A zero value stays 0 in every sample (clipped mean 0).
        for s in a { #expect(s.values[5] == 0) }
        for bad in [0, -1, Int.max, Int.max / 6] { #expect(throws: SimulatoError.self) { _ = try m.poissonNoise(dose: 3, samples: bad, seed: 1) } }
        #expect(throws: SimulatoError.self) { _ = try m.poissonNoise(dose: -1, samples: 2, seed: 1) }
        let perArea = try m.poissonNoise(dosePerArea: 1, samplingX: 0.5, samplingY: 4, samples: 2, seed: 9)
        let direct = try m.poissonNoise(dose: 2, samples: 2, seed: 9)
        #expect(perArea.map(\.values) == direct.map(\.values))
        #expect(throws: SimulatoError.self) { _ = try m.poissonNoise(dosePerArea: 1, samplingX: 0, samplingY: 4, samples: 2, seed: 9) }
    }

    @Test func frozenPhononWavesKeepOneScanPerConfiguration() throws {
        let (probe, potential) = try WavesDetectorRefusalTests.setup()
        let scan = try WavesDetectorRefusalTests.scan()
        let structure = try aluminiumStructure(cells: (1, 1, 2))
        let phonons = try FrozenPhonons(structure: structure, sigmas: [13: 0.1], configurations: 2, seed: 5)
        let potentials = try phonons.lazyPotentials(nx: 48, ny: 48, sliceThickness: 2.025,
                                                    parameterization: LobatoParameterization(piConstant: .abTEMFloat32))
        let detector = try WavesDetector(rows: 16, columns: 16)
        let ensemble = try STEM.scanFrozenPhonons(probe: probe, potentials: potentials, scan: scan, detector: detector)
        #expect(ensemble.count == 2)
        // Each configuration is the single-potential scan of that configuration, bit for bit.
        for (k, p) in potentials.enumerated() {
            let single = try STEM.scan(probe: probe, potentials: [p], scan: scan, detector: detector)
            #expect(ensemble[k].scanShape == single.scanShape)
            #expect(ensemble[k].real == single.real && ensemble[k].imaginary == single.imaginary)
        }
        // The configurations differ (displaced atoms), so nothing was averaged.
        #expect(ensemble[0].real != ensemble[1].real)
        // Built potentials give the same waves.
        let built = try STEM.scanFrozenPhonons(probe: probe, potentials: potentials.map { try $0.build() }, scan: scan,
                                               detector: detector)
        for k in 0..<2 { #expect(built[k].real == ensemble[k].real && built[k].imaginary == ensemble[k].imaginary) }
        // A configuration on another grid is refused before any work.
        let other = try Potential(structure: structure, nx: 40, ny: 40, sliceThickness: 2.025)
        #expect(throws: SimulatoError.gridMismatch) {
            try STEM.scanFrozenPhonons(probe: probe, potentials: [potentials[0], other], scan: scan, detector: detector)
        }
        // One static potential is an ensemble of one.
        #expect(try STEM.scanFrozenPhonons(probe: probe, potentials: [potential], scan: scan, detector: detector).count == 1)
    }

    @Test func frozenPhononWavesCountEveryConfigurationAgainstTheMemoryLimit() throws {
        let (probe, potential) = try WavesDetectorRefusalTests.setup()
        let scan = try WavesDetectorRefusalTests.scan()
        let detector = try WavesDetector(rows: 10, columns: 11)
        let one = WavesDetector.requiredBytes(positions: 6, rows: 10, columns: 11)
        #expect(throws: WavesDetectorError.tooLarge(positions: 18, rows: 10, columns: 11, required: 3 * one, limit: 3 * one - 1)) {
            try STEM.scanFrozenPhonons(probe: probe, potentials: [potential, potential, potential], scan: scan,
                                       detector: detector, memoryLimit: 3 * one - 1)
        }
        #expect(try STEM.scanFrozenPhonons(probe: probe, potentials: [potential, potential, potential], scan: scan,
                                           detector: detector, memoryLimit: 3 * one).count == 3)
        #expect(throws: SimulatoError.self) {
            try STEM.scanFrozenPhonons(probe: probe, potentials: [Potential](), scan: scan, detector: detector)
        }
    }
}
