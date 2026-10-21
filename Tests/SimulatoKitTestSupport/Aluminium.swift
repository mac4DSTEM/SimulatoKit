//
//  Aluminium.swift
//  Role: the face-centred aluminium cell, repeated, that the unit tests and the parity tests both build.
//

import SimulatoCore

package func aluminiumStructure(cells: (Int, Int, Int)) throws -> AtomicStructure {
    let cell = try AtomicStructure(box: SIMD3(4.05, 4.05, 4.05), atoms: [
        .init(atomicNumber: 13, position: SIMD3(0, 0, 0)), .init(atomicNumber: 13, position: SIMD3(0, 2.025, 2.025)),
        .init(atomicNumber: 13, position: SIMD3(2.025, 0, 2.025)), .init(atomicNumber: 13, position: SIMD3(2.025, 2.025, 0)),
    ])
    return cell.repeated(cells.0, cells.1, cells.2)
}
