//
//  Slicing.swift
//  Role: the slices of a structure along z (registration section 3.4; abTEM `_validate_slice_thickness`,
//        `SliceIndexedAtoms`, `_prepare_atoms`): periodic wrap, the top-face rule, equal slices from a requested
//        thickness, and the nudged cumulative edges.
//

import Foundation
import SimulatoCore

/// How a structure is cut into slices along the beam.
public enum Slicing {
    /// The most slices one potential may have (a guard against a requested thickness so small that the count overflows).
    public static let maximumSliceCount = 1_000_000

    /// n = ⌈Lz / s⌉ slices of Lz / n each (abTEM's rule for a single requested thickness).
    public static func thicknesses(depth: Double, requested s: Double) throws -> [Double] {
        guard s.isFinite, s > 0 else { throw SimulatoError.invalidSliceThickness(s) }
        let n = (depth / s).rounded(.up)
        // Refuter finding (2026-10-15): a tiny s overflowed Int(n) and trapped. More than a million slices is refused.
        guard n.isFinite, n >= 1, n <= Double(Slicing.maximumSliceCount) else { throw SimulatoError.invalidSliceThickness(s) }
        return [Double](repeating: depth / n, count: Int(n))
    }

    /// A list of thicknesses, which must sum to the depth within NumPy's `isclose` (rtol 1e−5, atol 1e−8).
    public static func validated(_ list: [Double], depth: Double) throws -> [Double] {
        guard !list.isEmpty, list.count <= maximumSliceCount else { throw SimulatoError.invalidSliceThickness(0) }
        for s in list where !(s.isFinite && s > 0) { throw SimulatoError.invalidSliceThickness(s) }
        let sum = list.reduce(0, +)
        guard abs(sum - depth) <= 1e-8 + 1e-5 * abs(depth) else {
            throw SimulatoError.sliceThicknessesDoNotSumToBox(sum: sum, depth: depth)
        }
        return list
    }

    /// The atoms wrapped into the box (fractional coordinate mod 1, as ASE's `wrap(eps=0)`), with z set to 0 for an atom
    /// within 1e−10 Å of the top face.
    public static func wrapped(_ structure: AtomicStructure) -> [AtomicStructure.Atom] {
        let box = structure.box
        return structure.atoms.map { atom in
            var p = atom.position
            for d in 0..<3 {
                var f = (p[d] / box[d]).truncatingRemainder(dividingBy: 1.0)
                if f < 0 { f += 1.0 }
                if f == 1.0 { f = 0 }
                p[d] = f * box[d]
            }
            if p.z > box.z - 1e-10 { p.z = 0 }
            return AtomicStructure.Atom(atomicNumber: atom.atomicNumber, position: p)
        }
    }

    /// The slice index of each z: e_{j−1} ≤ z < e_j with e_j = Σ_{i ≤ j} s_i − 1e−12 (NumPy's `digitize`); a z past the
    /// last edge goes into the last slice.
    public static func sliceIndices(z: [Double], thicknesses: [Double]) -> [Int] {
        var edges: [Double] = []
        var acc = 0.0
        for s in thicknesses {
            acc += s
            edges.append(acc)
        }
        edges = edges.map { $0 - 1e-12 }
        return z.map { value in
            var lo = 0, hi = edges.count
            while lo < hi {  // first edge strictly greater than value
                let mid = (lo + hi) / 2
                if edges[mid] <= value { lo = mid + 1 } else { hi = mid }
            }
            return min(lo, thicknesses.count - 1)
        }
    }
}
