//
//  AtomicStructure.swift
//  Role: the atoms to simulate: an orthogonal box and atoms with Cartesian positions in Å.
//        MARKED FOR REPLACEMENT: a minimal local type until DiffractoKit has its first tag; then SimulatoCore takes
//        DiffractoKit's `CrystalStructure` (cell, sites, symmetry) and keeps only the orthogonal-box view multislice
//        needs (registration section 1).
//

/// An orthogonal box with atoms in Cartesian coordinates [Å], the beam along z.
public struct AtomicStructure: Sendable, Hashable {
    public struct Atom: Sendable, Hashable {
        /// Atomic number Z.
        public var atomicNumber: Int
        /// Cartesian position (x, y, z) [Å].
        public var position: SIMD3<Double>

        public init(atomicNumber: Int, position: SIMD3<Double>) {
            self.atomicNumber = atomicNumber
            self.position = position
        }
    }

    /// The box's sides (Lx, Ly, Lz) [Å]; the cell is diagonal.
    public let box: SIMD3<Double>
    public let atoms: [Atom]

    /// - Throws: `SimulatoError.invalidBox` for a non-finite or non-positive side, `invalidPosition` for a non-finite
    ///   coordinate. Atoms outside the box are wrapped into it when the potential is built (periodic boundaries).
    public init(box: SIMD3<Double>, atoms: [Atom]) throws {
        guard box.x.isFinite, box.y.isFinite, box.z.isFinite, box.x > 0, box.y > 0, box.z > 0 else {
            throw SimulatoError.invalidBox(box.x, box.y, box.z)
        }
        for (i, a) in atoms.enumerated()
        where !(a.position.x.isFinite && a.position.y.isFinite && a.position.z.isFinite) {
            throw SimulatoError.invalidPosition(index: i)
        }
        self.box = box
        self.atoms = atoms
    }

    /// The structure repeated nx × ny × nz times, in ASE's order (x outermost, then y, then z, the atoms innermost).
    public func repeated(_ nx: Int, _ ny: Int, _ nz: Int) -> AtomicStructure {
        precondition(nx >= 1 && ny >= 1 && nz >= 1, "repeat counts must be positive")
        var out: [Atom] = []
        out.reserveCapacity(atoms.count * nx * ny * nz)
        for mx in 0..<nx {
            for my in 0..<ny {
                for mz in 0..<nz {
                    let shift = SIMD3<Double>(Double(mx) * box.x, Double(my) * box.y, Double(mz) * box.z)
                    for a in atoms { out.append(Atom(atomicNumber: a.atomicNumber, position: a.position + shift)) }
                }
            }
        }
        // The box sides are finite and positive, and the shifted positions finite: the init cannot throw.
        return try! AtomicStructure(box: box * SIMD3(Double(nx), Double(ny), Double(nz)), atoms: out)
    }
}
