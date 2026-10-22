//
//  IntegratedCenterOfMass.swift
//  Role: the centre-of-mass image of a scan and its integral (abTEM `DiffractionPatterns.integrated_center_of_mass`, i.e.
//        `center_of_mass().integrate_gradient()`): the two values of a `CenterOfMassDetector` over the scan are the x and
//        y components of a gradient field, integrated in Fourier space. abTEM applies no further scaling: the image is in
//        the units of the centre of mass times Å, shifted to a minimum of 0 (its label is "iCOM", units "arb. unit").
//

import SimulatoCore

extension ScanMeasurement {
    /// The x and y components of a two-value measurement (a `CenterOfMassDetector`'s) as two images over the scan, with the
    /// scan's pixel sizes `samplingX` and `samplingY` [Å].
    /// - Throws: `SimulatoError.invalidScan` unless every position holds exactly two values and the scan has two axes.
    public func gradientImages(samplingX: Double, samplingY: Double) throws -> (x: RealField, y: RealField) {
        guard valueShape == [2] else {
            throw SimulatoError.invalidScan("a centre-of-mass image needs two values per position (x, y), this measurement has shape \(valueShape)")
        }
        guard scanShape.count == 2 else {
            throw SimulatoError.invalidScan("a centre-of-mass image needs a two-axis scan, this measurement has scan shape \(scanShape)")
        }
        let n = scanShape[0] * scanShape[1]
        var gx = [Double](repeating: 0, count: n), gy = gx
        for i in 0..<n {
            gx[i] = values[2 * i]
            gy[i] = values[2 * i + 1]
        }
        return (RealField(nx: scanShape[0], ny: scanShape[1], values: gx, samplingX: samplingX, samplingY: samplingY),
                RealField(nx: scanShape[0], ny: scanShape[1], values: gy, samplingX: samplingX, samplingY: samplingY))
    }

    /// The integrated centre of mass (iCOM; the integrated DPC of a two-value DPC measurement is the same operation): the
    /// Fourier-space integral of the centre-of-mass image (`RealField.integratingGradient`), in the
    /// units of the centre of mass times Å, minimum 0, sampled as the scan.
    /// - Throws: `SimulatoError.invalidScan` unless every position holds two values, `invalidSampling` for a
    ///   non-finite or non-positive sampling, `unsupportedGridLength` for a scan axis shorter than 2.
    public func integratedCenterOfMass(samplingX: Double, samplingY: Double) throws -> RealField {
        let g = try gradientImages(samplingX: samplingX, samplingY: samplingY)
        return try RealField.integratingGradient(x: g.x, y: g.y)
    }
}
