//
//  SimulatoError.swift
//  Role: every refusal the kit makes, named; messages say what was wrong and what is accepted.
//

/// A refusal by SimulatoKit.
public enum SimulatoError: Error, Sendable, Equatable, CustomStringConvertible {
    case invalidEnergy(Double)
    case unsupportedGridLength(Int)
    case invalidExtent(Double, Double)
    case invalidSampling(Double, Double)
    case invalidTilt(String)
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
    case invalidFrozenPhonons(String)
    case invalidScan(String)
    case invalidDetector(String)
    /// A resampled wave of `rows` × `columns` points needs `required` bytes, more than the `limit` allowed (`Int.max`
    /// bytes for a size that overflows).
    case gridTooLarge(rows: Int, columns: Int, required: Int, limit: Int)

    public var description: String {
        switch self {
        case .invalidEnergy(let e): return "electron energy must be finite and positive, got \(e) eV"
        case .unsupportedGridLength(let n):
            return "grid length \(n) is not supported: lengths must be at least 2"
        case .invalidExtent(let x, let y): return "grid extent must be finite and positive, got (\(x), \(y)) Å"
        case .invalidSampling(let x, let y): return "grid sampling must be finite and positive and give fewer than 1e9 points, got (\(x), \(y)) Å"
        case .invalidTilt(let why): return "invalid tilt: \(why)"
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
        case .invalidFrozenPhonons(let why): return "invalid frozen phonons: \(why)"
        case .invalidScan(let why): return "invalid scan: \(why)"
        case .invalidDetector(let why): return "invalid detector: \(why)"
        case .gridTooLarge(let r, let c, let need, let limit):
            return "a wave of \(r) x \(c) points needs \(need) bytes, more than the \(limit) bytes allowed (half of the physical memory)"
        }
    }
}
