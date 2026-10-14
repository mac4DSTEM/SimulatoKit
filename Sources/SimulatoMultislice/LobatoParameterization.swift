//
//  LobatoParameterization.swift
//  Role: the Lobato & Van Dyck (2014) parameterisation, abTEM's default (registration section 3.3): the scaled
//        parameters A = π² a / b^{3/2} / κ, B = 2π / √b, and the projected scattering factor of an infinite
//        projection, f(k²) = 8π Σ [A/B/(q + B²) + A B/(q + B²)²], q = 4π² k², summed in abTEM's order.
//

import Foundation
import SimulatoCore

/// The Lobato & Van Dyck electron scattering parameterisation (Acta Cryst. A 70, 636–649, 2014), Z = 1 ... 103.
public struct LobatoParameterization: Sendable, Hashable {
    /// The π in f(k²). abTEM 1.0.10 evaluates it as a float32 constant, also in float64 mode, which moves f by about
    /// 5e−8 relative; `.exact` is the default, `.abTEMFloat32` reproduces abTEM's numbers (registration 3.3).
    public enum PiConstant: Sendable, Hashable {
        case exact
        case abTEMFloat32
    }

    public var piConstant: PiConstant

    public init(piConstant: PiConstant = .exact) {
        self.piConstant = piConstant
    }

    /// Whether Z has parameters.
    public static func supports(atomicNumber z: Int) -> Bool {
        z >= 1 && z <= lobatoTable.count && lobatoTable[z - 1] != nil
    }

    /// The table's (a, b) for Z, as abTEM ships them.
    public static func parameters(atomicNumber z: Int) throws -> (a: [Double], b: [Double]) {
        guard supports(atomicNumber: z), let p = lobatoTable[z - 1] else { throw SimulatoError.unsupportedElement(z) }
        return p
    }

    /// The scaled parameters (A, B) of the projected scattering factor (abTEM's `scaled_parameters`, exact π).
    public static func scaledParameters(atomicNumber z: Int) throws -> (A: [Double], B: [Double]) {
        let p = try parameters(atomicNumber: z)
        let pi2 = Double.pi * Double.pi
        let a = (0..<5).map { i in pi2 * p.a[i] / pow(p.b[i], 1.5) / PhysicalConstants.kappa }
        let b = (0..<5).map { i in 2 * Double.pi / p.b[i].squareRoot() }
        return (a, b)
    }

    /// The float32-rounded π and π² abTEM uses in f (3.1415927410125732, 9.869604110717773).
    public static let float32Pi = Double(Float(Double.pi))
    public static let float32PiSquared = Double(Float(Double.pi * Double.pi))

    /// f(k²) [V·Å²] for one element at the squared spatial frequencies k² [1/Å²].
    public func projectedScatteringFactor(atomicNumber z: Int, k2: [Double]) throws -> [Double] {
        let s = try Self.scaledParameters(atomicNumber: z)
        let pi: Double, pi2: Double
        switch piConstant {
        case .exact: (pi, pi2) = (Double.pi, Double.pi * Double.pi)
        case .abTEMFloat32: (pi, pi2) = (Self.float32Pi, Self.float32PiSquared)
        }
        let fourPi2 = 4 * pi2
        let eightPi = 8 * pi
        return k2.map { k in
            let q = fourPi2 * k
            var sum = 0.0
            for i in 0..<5 {
                let d = q + s.B[i] * s.B[i]
                sum = sum + (s.A[i] / s.B[i] / d + s.A[i] * s.B[i] / (d * d))
            }
            return eightPi * sum
        }
    }

    /// f(k²) at one k².
    public func projectedScatteringFactor(atomicNumber z: Int, k2: Double) throws -> Double {
        try projectedScatteringFactor(atomicNumber: z, k2: [k2])[0]
    }
}
