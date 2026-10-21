//
//  Comparison.swift
//  Role: the relative-difference helpers both the unit tests and the parity tests use (nothing is copied).
//

import SimulatoCore

/// max |a − b| / max |b| (0 when both are zero).
package func relMax(_ a: [Double], _ b: [Double]) -> Double {
    precondition(a.count == b.count, "relMax: \(a.count) vs \(b.count) values")
    var d = 0.0, m = 0.0
    for i in a.indices {
        d = max(d, abs(a[i] - b[i]))
        m = max(m, abs(b[i]))
    }
    return m > 0 ? d / m : d
}

/// max |a − b| / max |b| for complex fields.
package func relMax(_ a: ComplexField, _ b: ComplexField) -> Double {
    precondition(a.nx == b.nx && a.ny == b.ny)
    var d = 0.0, m = 0.0
    for i in a.real.indices {
        let dr = a.real[i] - b.real[i], di = a.imaginary[i] - b.imaginary[i]
        d = max(d, (dr * dr + di * di).squareRoot())
        m = max(m, (b.real[i] * b.real[i] + b.imaginary[i] * b.imaginary[i]).squareRoot())
    }
    return m > 0 ? d / m : d
}

package func relDiff(_ a: Double, _ b: Double) -> Double { b == 0 ? abs(a) : abs(a - b) / abs(b) }
