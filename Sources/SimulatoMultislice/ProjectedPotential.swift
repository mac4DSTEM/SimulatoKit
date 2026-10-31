//
//  ProjectedPotential.swift
//  Role: the projected potential of each slice, abTEM's default `Potential(parametrization="lobato",
//        projection="infinite")` (`integrals.py` `ScatteringFactorProjectionIntegrals`):
//        bilinear deltas per element, FFT, times f(k²) / S(k), inverse FFT, real part, summed over elements.
//

import Foundation
import SimulatoCore

/// One slice: its thickness [Å] and its projected potential V [V·Å] on the grid.
public struct PotentialSlice: Sendable, Hashable {
    public let thickness: Double
    public let values: [Double]

    /// A slice of thickness [Å] with V [V·Å] row-major on its grid (`ProjectedPotential(grid:slices:)` checks the count).
    public init(thickness: Double, values: [Double]) {
        self.thickness = thickness
        self.values = values
    }
}

/// The projected potential of a structure, slice by slice, on a grid.
public struct ProjectedPotential: Sendable, Hashable {
    public let grid: Grid
    public let slices: [PotentialSlice]
    /// The slice each atom of the structure was put in (in the structure's atom order).
    public let sliceOfAtom: [Int]

    /// Builds the potential with slices of `sliceThickness` (abTEM's rule: ⌈Lz / s⌉ equal slices), at most
    /// `maximumSliceCount`.
    public init(structure: AtomicStructure, nx: Int, ny: Int, sliceThickness: Double,
                parameterization: LobatoParameterization = LobatoParameterization(),
                maximumSliceCount: Int = Slicing.maximumSliceCount) throws {
        let t = try Slicing.thicknesses(depth: structure.box.z, requested: sliceThickness, maximumCount: maximumSliceCount)
        try self.init(structure: structure, nx: nx, ny: ny, sliceThicknesses: t, parameterization: parameterization,
                      maximumSliceCount: maximumSliceCount)
    }

    /// Builds the potential with the given slice thicknesses (they must sum to the box's depth): `Potential(...).build()`.
    public init(structure: AtomicStructure, nx: Int, ny: Int, sliceThicknesses: [Double],
                parameterization: LobatoParameterization = LobatoParameterization(),
                maximumSliceCount: Int = Slicing.maximumSliceCount) throws {
        self = try Potential(structure: structure, nx: nx, ny: ny, sliceThicknesses: sliceThicknesses,
                             parameterization: parameterization, maximumSliceCount: maximumSliceCount).build()
    }

    init(grid: Grid, slices: [PotentialSlice], sliceOfAtom: [Int]) {
        self.grid = grid
        self.slices = slices
        self.sliceOfAtom = sliceOfAtom
    }

    /// A potential from slices computed elsewhere (for example the slices of another potential).
    public init(grid: Grid, slices: [PotentialSlice]) {
        precondition(slices.allSatisfy { $0.values.count == grid.count }, "every slice must match the grid")
        self.grid = grid
        self.slices = slices
        self.sliceOfAtom = []
    }

    /// The total thickness [Å].
    public var thickness: Double { slices.reduce(0) { $0 + $1.thickness } }

    /// The sum of the slices as an image [V·Å] on the potential's grid (abTEM `PotentialArray.project`), added slice by
    /// slice in order.
    public func project() -> RealField {
        var sum = [Double](repeating: 0, count: grid.count)
        for slice in slices {
            for i in sum.indices { sum[i] += slice.values[i] }
        }
        return RealField(nx: grid.nx, ny: grid.ny, values: sum, samplingX: grid.samplingX, samplingY: grid.samplingY)
    }

    /// The potential tiled `repetitionsX` × `repetitionsY` times in the plane and `repetitionsZ` times along the beam
    /// (abTEM `PotentialArray.tile`): every slice tiled as an image, the extent multiplied and the sampling kept, then the
    /// whole stack of slices repeated `repetitionsZ` times. The result carries no `sliceOfAtom`, as abTEM's array
    /// carries no atoms.
    /// - Throws: `ImageOperationError.invalidRepetitions` (which reports the x and y counts only) for a count below 1 on
    ///   any axis, a tiled slice of more than 2^28 values, more than `Slicing.maximumSliceCount` slices or more than 2^31
    ///   values in all; all checked before anything is allocated.
    public func tiled(repetitionsX: Int, repetitionsY: Int, repetitionsZ: Int = 1) throws -> ProjectedPotential {
        let refusal = ImageOperationError.invalidRepetitions(repetitionsX, repetitionsY)
        guard repetitionsX >= 1, repetitionsY >= 1, repetitionsZ >= 1 else { throw refusal }
        let (sliceCount, o1) = slices.count.multipliedReportingOverflow(by: repetitionsZ)
        let (mx, o2) = grid.nx.multipliedReportingOverflow(by: repetitionsX)
        let (my, o3) = grid.ny.multipliedReportingOverflow(by: repetitionsY)
        let (perSlice, o4) = mx.multipliedReportingOverflow(by: my)
        let (total, o5) = perSlice.multipliedReportingOverflow(by: sliceCount)
        guard !o1, !o2, !o3, !o4, !o5, sliceCount <= Slicing.maximumSliceCount, perSlice <= 1 << 28, total <= 1 << 31 else {
            throw refusal
        }
        var tiledSlices: [PotentialSlice] = []
        tiledSlices.reserveCapacity(slices.count)
        for slice in slices {
            let image = RealField(nx: grid.nx, ny: grid.ny, values: slice.values,
                                  samplingX: grid.samplingX, samplingY: grid.samplingY)
            let tiledImage = try image.tiled(repetitionsX: repetitionsX, repetitionsY: repetitionsY)
            tiledSlices.append(PotentialSlice(thickness: slice.thickness, values: tiledImage.values))
        }
        let tiledGrid = try Grid(nx: mx, ny: my,
                                 extentX: grid.extentX * Double(repetitionsX), extentY: grid.extentY * Double(repetitionsY))
        let stack = Array([[PotentialSlice]](repeating: tiledSlices, count: repetitionsZ).joined())
        return ProjectedPotential(grid: tiledGrid, slices: stack, sliceOfAtom: [])
    }
}

/// The per-grid work of the projection: f(k²) / S(k) per element (cached) and the FFT plan.
final class SliceBuilder {
    let grid: Grid
    let parameterization: LobatoParameterization
    private let fft: FFT2D
    private let k2: [Double]
    private let compensation: [Double]
    private var kernels: [Int: [Double]] = [:]

    init(grid: Grid, parameterization: LobatoParameterization) throws {
        self.grid = grid
        self.parameterization = parameterization
        fft = try FFT2D(nx: grid.nx, ny: grid.ny)
        let kx = grid.frequenciesX, ky = grid.frequenciesY
        var k2 = [Double](repeating: 0, count: grid.count)
        for ix in 0..<grid.nx { for iy in 0..<grid.ny { k2[ix * grid.ny + iy] = kx[ix] * kx[ix] + ky[iy] * ky[iy] } }
        self.k2 = k2
        compensation = Self.compensation(grid: grid)
    }

    /// abTEM's `sinc` (integrals.py): S = sin(κ)/κ · Δx Δy with κ = √((kx Δx)² + (ky Δy)²) in cycles per pixel (no π),
    /// S(0) = Δx Δy.
    static func compensation(grid: Grid) -> [Double] {
        let kx = grid.frequenciesX, ky = grid.frequenciesY
        let dx = grid.samplingX, dy = grid.samplingY
        let dk2 = dx * dy
        var s = [Double](repeating: 0, count: grid.count)
        for ix in 0..<grid.nx {
            for iy in 0..<grid.ny {
                let a = kx[ix] * dx, b = ky[iy] * dy
                let k = (a * a + b * b).squareRoot()
                s[ix * grid.ny + iy] = (ix == 0 && iy == 0) ? dk2 : sin(k) / k * dk2
            }
        }
        return s
    }

    /// f(k²) / S(k) for element Z.
    func kernel(_ z: Int) throws -> [Double] {
        if let k = kernels[z] { return k }
        let f = try parameterization.projectedScatteringFactor(atomicNumber: z, k2: k2)
        var out = f
        for i in out.indices { out[i] = f[i] / compensation[i] }
        kernels[z] = out
        return out
    }

    /// V of one slice from its (wrapped) atoms: elements in increasing Z, each Re(IFFT(FFT(deltas) · f / S)).
    func potential(of atoms: [AtomicStructure.Atom]) throws -> [Double] {
        var v = [Double](repeating: 0, count: grid.count)
        let numbers = Set(atoms.map(\.atomicNumber)).sorted()
        for z in numbers {
            let positions = atoms.filter { $0.atomicNumber == z }.map { SIMD2($0.position.x / grid.samplingX,
                                                                                $0.position.y / grid.samplingY) }
            var field = ComplexField(nx: grid.nx, ny: grid.ny)
            Self.superposeDeltas(positions, into: &field.real, nx: grid.nx, ny: grid.ny)
            fft.forward(&field)
            let k = try kernel(z)
            for i in k.indices {
                field.real[i] *= k[i]
                field.imaginary[i] *= k[i]
            }
            fft.inverse(&field)
            for i in v.indices { v[i] += field.real[i] }
        }
        return v
    }

    /// abTEM's `superpose_deltas` (subpixel): weights 1 + xy − y − x, x − xy, y − xy, xy at (i, j), (i+1, j), (i, j+1),
    /// (i+1, j+1), indices mod the grid; accumulated weight by weight, atom by atom (NumPy's `add.at` order).
    static func superposeDeltas(_ positions: [SIMD2<Double>], into array: inout [Double], nx: Int, ny: Int) {
        let n = positions.count
        var rows = [Int](repeating: 0, count: n), cols = rows
        var w = [[Double]](repeating: [Double](repeating: 0, count: n), count: 4)
        for (a, p) in positions.enumerated() {
            let r = Int(p.x.rounded(.down)), c = Int(p.y.rounded(.down))
            let x = p.x - Double(r), y = p.y - Double(c)
            let xy = x * y
            rows[a] = r
            cols[a] = c
            w[0][a] = 1 + xy - y - x
            w[1][a] = x - xy
            w[2][a] = y - xy
            w[3][a] = xy
        }
        func mod(_ a: Int, _ m: Int) -> Int { let r = a % m; return r < 0 ? r + m : r }
        for k in 0..<4 {
            for a in 0..<n {
                let i = mod(rows[a] + (k == 1 || k == 3 ? 1 : 0), nx)
                let j = mod(cols[a] + (k >= 2 ? 1 : 0), ny)
                array[i * ny + j] += w[k][a]
            }
        }
    }
}
