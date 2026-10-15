//
//  Potential.swift
//  Role: abTEM's lazy `Potential` (`potentials/iam.py`; registration-performance.md 1.1): a structure cut into slices on
//        a grid, the wrapped atoms with their slice and the parameterisation, without the arrays. `build()` gives the
//        `ProjectedPotential` (abTEM's `PotentialArray`) on the CPU; on `.metal` the multislice builds each slice on the
//        GPU instead.
//

import SimulatoCore

/// A structure sliced on a grid, ready to be built into a projected potential (on the CPU by `build()`, or slice by slice
/// on the GPU by the `.metal` backend).
public struct Potential: Sendable, Hashable {
    public let grid: Grid
    /// The thickness of every slice [Å].
    public let sliceThicknesses: [Double]
    public let parameterization: LobatoParameterization
    /// The atoms wrapped into the box (registration-slice-1.md 3.4), in the structure's atom order.
    public let atoms: [AtomicStructure.Atom]
    /// The slice of each atom, in the same order.
    public let sliceOfAtom: [Int]

    /// Slices of `sliceThickness` (abTEM's rule: ⌈Lz / s⌉ equal slices).
    public init(structure: AtomicStructure, nx: Int, ny: Int, sliceThickness: Double,
                parameterization: LobatoParameterization = LobatoParameterization()) throws {
        let t = try Slicing.thicknesses(depth: structure.box.z, requested: sliceThickness)
        try self.init(structure: structure, nx: nx, ny: ny, sliceThicknesses: t, parameterization: parameterization)
    }

    /// The given slice thicknesses (they must sum to the box's depth).
    public init(structure: AtomicStructure, nx: Int, ny: Int, sliceThicknesses: [Double],
                parameterization: LobatoParameterization = LobatoParameterization()) throws {
        let grid = try Grid(nx: nx, ny: ny, extentX: structure.box.x, extentY: structure.box.y)
        let thicknesses = try Slicing.validated(sliceThicknesses, depth: structure.box.z)
        for a in structure.atoms where !LobatoParameterization.supports(atomicNumber: a.atomicNumber) {
            throw SimulatoError.unsupportedElement(a.atomicNumber)
        }
        let atoms = Slicing.wrapped(structure)
        self.grid = grid
        self.sliceThicknesses = thicknesses
        self.parameterization = parameterization
        self.atoms = atoms
        self.sliceOfAtom = Slicing.sliceIndices(z: atoms.map(\.position.z), thicknesses: thicknesses)
    }

    /// The number of slices.
    public var sliceCount: Int { sliceThicknesses.count }

    /// The total thickness [Å].
    public var thickness: Double { sliceThicknesses.reduce(0, +) }

    /// The atoms of slice j, in the structure's order.
    public func atoms(inSlice j: Int) -> [AtomicStructure.Atom] {
        atoms.indices.filter { sliceOfAtom[$0] == j }.map { atoms[$0] }
    }

    /// The projected potential of every slice, on the CPU in float64 (registration-slice-1.md 3.5).
    public func build() throws -> ProjectedPotential {
        let builder = try SliceBuilder(grid: grid, parameterization: parameterization)
        var slices: [PotentialSlice] = []
        slices.reserveCapacity(sliceThicknesses.count)
        for (j, dz) in sliceThicknesses.enumerated() {
            slices.append(PotentialSlice(thickness: dz, values: try builder.potential(of: atoms(inSlice: j))))
        }
        return ProjectedPotential(grid: grid, slices: slices, sliceOfAtom: sliceOfAtom)
    }
}
