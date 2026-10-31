//
//  ScanOperations.swift
//  Role: what abTEM does to a finished scan: the Gaussian source size (`gaussian_source_size`: scipy's `gaussian_filter`
//        over the two scan axes, mode "wrap") and Poisson noise (`poisson_noise`: values times the dose, clipped at 0,
//        one Poisson draw per value; `samples` draws as an array).
//

import Foundation
import SimulatoCore
import SimulatoMultislice

extension ScanMeasurement {
    /// The scan with the finite source size of a Gaussian probe of standard deviation `sigma` [Å] (abTEM
    /// `gaussian_source_size`): every value filtered over the scan, with the scan's pixel size `samplingX` and
    /// `samplingY` [Å] (the scan's extent over its points). The kernel is scipy's: radius ⌊4 σ / d + ½⌋ pixels,
    /// normalised, periodic. The intensity is mixed across scan positions, so integrating the patterns before or after
    /// gives the same image.
    /// - Throws: `invalidScan` for a negative or non-finite σ or sampling.
    public func gaussianSourceSize(sigmaX: Double, sigmaY: Double, samplingX: Double, samplingY: Double) throws -> ScanMeasurement {
        guard sigmaX.isFinite, sigmaY.isFinite, sigmaX >= 0, sigmaY >= 0 else {
            throw SimulatoError.invalidScan("the source size must be finite and non-negative, got (\(sigmaX), \(sigmaY)) Å")
        }
        guard samplingX.isFinite, samplingY.isFinite, samplingX > 0, samplingY > 0 else {
            throw SimulatoError.invalidScan("the scan sampling must be finite and positive, got (\(samplingX), \(samplingY)) Å")
        }
        guard scanShape.count == 2 else { throw SimulatoError.invalidScan("the source size filters a two-axis scan; a line scan has one axis") }
        let nx = scanShape[0], ny = scanShape[1], v = valueCount
        var data = values
        // Along y (the second scan axis) then x, each with its own kernel; scipy filters the axes in order.
        Self.filterAxis(&data, count: nx * ny * v, length: ny, stride: v, outer: nx, outerStride: ny * v, inner: v,
                        sigma: sigmaY / samplingY)
        Self.filterAxis(&data, count: nx * ny * v, length: nx, stride: ny * v, outer: 1, outerStride: 0, inner: ny * v,
                        sigma: sigmaX / samplingX)
        return ScanMeasurement(scanShape: scanShape, valueShape: valueShape, values: data)
    }

    /// The same with one σ for both axes.
    public func gaussianSourceSize(sigma: Double, samplingX: Double, samplingY: Double) throws -> ScanMeasurement {
        try gaussianSourceSize(sigmaX: sigma, sigmaY: sigma, samplingX: samplingX, samplingY: samplingY)
    }

    /// scipy's `gaussian_filter1d` weights: exp(−x²/2σ²) for |x| ≤ ⌊4σ + ½⌋, summing to 1.
    static func gaussianKernel(sigmaInPixels s: Double) -> [Double] {
        let radius = Int(4.0 * s + 0.5)
        var w = (-radius...radius).map { exp(-0.5 * Double($0 * $0) / (s * s)) }
        let sum = w.reduce(0, +)
        for i in w.indices { w[i] /= sum }
        return w
    }

    private static func filterAxis(_ data: inout [Double], count: Int, length: Int, stride: Int, outer: Int, outerStride: Int,
                                   inner: Int, sigma: Double) {
        guard sigma > 0 else { return }
        let kernel = gaussianKernel(sigmaInPixels: sigma)
        let radius = kernel.count / 2
        var line = [Double](repeating: 0, count: length)
        for o in 0..<outer {
            for i in 0..<inner {
                let base = o * outerStride + i
                for k in 0..<length { line[k] = data[base + k * stride] }
                for k in 0..<length {
                    var acc = 0.0
                    for (j, w) in kernel.enumerated() {
                        var idx = (k + j - radius) % length
                        if idx < 0 { idx += length }
                        acc += w * line[idx]
                    }
                    data[base + k * stride] = acc
                }
            }
        }
    }

    /// Poisson noise for a total `dose` [electrons per pattern] (abTEM `poisson_noise(total_dose=)`): every value times the
    /// dose, clipped at 0, then drawn from a Poisson distribution with that mean. The generator is the kit's own (a seed
    /// gives the same noise on every run, not the bits of NumPy's).
    /// - Throws: `invalidScan` for a negative or non-finite dose.
    public func poissonNoise(dose: Double, seed: UInt64) throws -> ScanMeasurement {
        guard dose.isFinite, dose >= 0 else {
            throw SimulatoError.invalidScan("the dose must be finite and non-negative, got \(dose)")
        }
        var rng = GaussianGenerator(seed: seed)
        return ScanMeasurement(scanShape: scanShape, valueShape: valueShape,
                               values: values.map { Double(PoissonSampler.draw(mean: max($0 * dose, 0), &rng)) })
    }

    /// Poisson noise for a dose per scan area [electrons per Å²] (abTEM `poisson_noise(dose_per_area=)`): the dose per
    /// pattern is the area of one scan pixel times `dosePerArea`.
    public func poissonNoise(dosePerArea: Double, samplingX: Double, samplingY: Double, seed: UInt64) throws -> ScanMeasurement {
        guard samplingX.isFinite, samplingY.isFinite, samplingX > 0, samplingY > 0 else {
            throw SimulatoError.invalidScan("the scan sampling must be finite and positive, got (\(samplingX), \(samplingY)) Å")
        }
        return try poissonNoise(dose: samplingX * samplingY * dosePerArea, seed: seed)
    }

    /// `samples` independent Poisson draws of the scan for a total `dose` [electrons per pattern] (abTEM
    /// `poisson_noise(total_dose=, samples=)`, whose sample axis is the array's index here): each sample is the scan
    /// times the dose, clipped at 0, drawn as `poissonNoise(dose:seed:)`. The samples continue one generator stream, so
    /// sample 0 is the single draw of the same seed.
    /// - Throws: `invalidScan` for fewer than one sample or a negative or non-finite dose.
    public func poissonNoise(dose: Double, samples: Int, seed: UInt64) throws -> [ScanMeasurement] {
        guard samples >= 1 else { throw SimulatoError.invalidScan("samples must be at least 1, got \(samples)") }
        // Refuse before any work when the samples would not fit in half of the physical memory (8 bytes a value).
        let (count, o1) = samples.multipliedReportingOverflow(by: values.count)
        let (bytes, o2) = count.multipliedReportingOverflow(by: 8)
        guard !o1, !o2, bytes <= Int(ProcessInfo.processInfo.physicalMemory / 2) else {
            throw SimulatoError.invalidScan("\(samples) samples of \(values.count) values exceed half of the physical memory")
        }
        guard dose.isFinite, dose >= 0 else {
            throw SimulatoError.invalidScan("the dose must be finite and non-negative, got \(dose)")
        }
        var rng = GaussianGenerator(seed: seed)
        let means = values.map { max($0 * dose, 0) }
        return (0..<samples).map { _ in
            ScanMeasurement(scanShape: scanShape, valueShape: valueShape,
                            values: means.map { Double(PoissonSampler.draw(mean: $0, &rng)) })
        }
    }

    /// The same for a dose per scan area [electrons per Å²] (abTEM `poisson_noise(dose_per_area=, samples=)`).
    public func poissonNoise(dosePerArea: Double, samplingX: Double, samplingY: Double, samples: Int,
                             seed: UInt64) throws -> [ScanMeasurement] {
        guard samplingX.isFinite, samplingY.isFinite, samplingX > 0, samplingY > 0 else {
            throw SimulatoError.invalidScan("the scan sampling must be finite and positive, got (\(samplingX), \(samplingY)) Å")
        }
        return try poissonNoise(dose: samplingX * samplingY * dosePerArea, samples: samples, seed: seed)
    }
}

/// Poisson draws: Knuth's product method below a mean of 10, Hörmann's PTRS (transformed rejection with squeeze) above.
enum PoissonSampler {
    static func uniform(_ rng: inout GaussianGenerator) -> Double {
        Double(rng.nextUInt64() >> 11) * (1.0 / 9007199254740992.0)
    }

    static func draw(mean lambda: Double, _ rng: inout GaussianGenerator) -> Int {
        if lambda <= 0 { return 0 }
        if lambda < 10 {
            let limit = exp(-lambda)
            var product = uniform(&rng), k = 0
            while product > limit {
                product *= uniform(&rng)
                k += 1
            }
            return k
        }
        let slam = lambda.squareRoot(), loglam = log(lambda)
        let b = 0.931 + 2.53 * slam
        let a = -0.059 + 0.02483 * b
        let invalpha = 1.1239 + 1.1328 / (b - 3.4)
        let vr = 0.9277 - 3.6224 / (b - 2)
        while true {
            let u = uniform(&rng) - 0.5
            let v = uniform(&rng)
            let us = 0.5 - abs(u)
            if us <= 0 { continue }
            let kd = ((2 * a / us + b) * u + lambda + 0.43).rounded(.down)
            guard kd.isFinite, abs(kd) < 9e18 else { continue }
            let k = Int(kd)
            if us >= 0.07 && v <= vr { return k }
            if k < 0 || (us < 0.013 && v > us) { continue }
            if log(v) + log(invalpha) - log(a / (us * us) + b) <= -lambda + Double(k) * loglam - lgamma(Double(k) + 1) {
                return k
            }
        }
    }
}
