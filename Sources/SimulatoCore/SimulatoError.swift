//
//  SimulatoError.swift
//  Role: every refusal the kit makes, named; messages say what was wrong and what is accepted.
//

/// A refusal by SimulatoKit.
public enum SimulatoError: Error, Sendable, Equatable, CustomStringConvertible {
    case invalidEnergy(Double)
    case unsupportedGridLength(Int)
    case invalidExtent(Double, Double)
    case invalidBox(Double, Double, Double)
    case unsupportedElement(Int)
    case invalidPosition(index: Int)
    case invalidSliceThickness(Double)
    case sliceThicknessesDoNotSumToBox(sum: Double, depth: Double)
    case invalidExitPlanes(String)
    case gridMismatch
    case invalidAperture(Double)
    case invalidAberration(String)
    case metalUnavailable(String)
    case metalFailure(String)

    public var description: String {
        switch self {
        case .invalidEnergy(let e): return "electron energy must be finite and positive, got \(e) eV"
        case .unsupportedGridLength(let n):
            return "grid length \(n) is not supported: lengths are f * 2^m with f in {1, 3, 5, 15} and m >= 3"
        case .invalidExtent(let x, let y): return "grid extent must be finite and positive, got (\(x), \(y)) Å"
        case .invalidBox(let x, let y, let z): return "the box must be orthogonal with finite positive sides, got (\(x), \(y), \(z)) Å"
        case .unsupportedElement(let z): return "element Z = \(z) has no Lobato parameters (Z = 1 ... 103 are tabulated)"
        case .invalidPosition(let i): return "atom \(i) has a non-finite position"
        case .invalidSliceThickness(let s): return "slice thickness must be finite and positive, got \(s) Å"
        case .sliceThicknessesDoNotSumToBox(let s, let d): return "slice thicknesses sum to \(s) Å, the box is \(d) Å deep"
        case .invalidExitPlanes(let why): return "invalid exit planes: \(why)"
        case .gridMismatch: return "the waves and the potential are on different grids"
        case .invalidAperture(let a): return "semiangle cutoff must be finite and positive (or infinite), got \(a) mrad"
        case .invalidAberration(let why): return "invalid aberration: \(why)"
        case .metalUnavailable(let why): return "the Metal backend is unavailable: \(why)"
        case .metalFailure(let why): return "the Metal backend failed: \(why)"
        }
    }
}
