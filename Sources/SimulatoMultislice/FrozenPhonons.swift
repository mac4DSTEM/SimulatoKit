//
//  FrozenPhonons.swift
//  Role: the thermal ensemble (abTEM `FrozenPhonons`): configurations of one structure, each atom displaced by a
//        Gaussian of standard deviation σ_Z along the chosen axes. Displacements are either given (e.g. abTEM's, from
//        the pins) or drawn by the kit's own seeded generator (SplitMix64 → xoshiro256**, Box–Muller), tested by its
//        statistics, not by abTEM's bits. The displaced structure is wrapped by the potential.
//

import Foundation
import SimulatoCore

/// A seeded standard-normal generator: SplitMix64 seeds xoshiro256**, Box–Muller turns pairs of uniforms into normals
/// (both outputs used). The same seed gives the same stream on every platform.
public struct GaussianGenerator: Sendable {
    private var s: (UInt64, UInt64, UInt64, UInt64)
    private var spare: Double?

    public init(seed: UInt64) {
        var x = seed
        func splitmix() -> UInt64 {
            x &+= 0x9E37_79B9_7F4A_7C15
            var z = x
            z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
            z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
            return z ^ (z >> 31)
        }
        s = (splitmix(), splitmix(), splitmix(), splitmix())
        spare = nil
    }

    /// xoshiro256** (Blackman and Vigna, public domain).
    public mutating func nextUInt64() -> UInt64 {
        func rotl(_ v: UInt64, _ k: UInt64) -> UInt64 { (v << k) | (v >> (64 - k)) }
        let result = rotl(s.1 &* 5, 7) &* 9
        let t = s.1 << 17
        s.2 ^= s.0
        s.3 ^= s.1
        s.1 ^= s.2
        s.0 ^= s.3
        s.2 ^= t
        s.3 = rotl(s.3, 45)
        return result
    }

    /// A standard normal deviate.
    public mutating func nextNormal() -> Double {
        if let z = spare {
            spare = nil
            return z
        }
        let u1 = (Double(nextUInt64() >> 11) + 1) * 0x1p-53  // (0, 1]
        let u2 = Double(nextUInt64() >> 11) * 0x1p-53  // [0, 1)
        let r = (-2 * log(u1)).squareRoot()
        let theta = 2 * Double.pi * u2
        spare = r * sin(theta)
        return r * cos(theta)
    }
}

/// A frozen-phonon ensemble of one structure.
public struct FrozenPhonons: Sendable, Hashable {
    /// The axes along which atoms are displaced.
    public struct Directions: OptionSet, Sendable, Hashable {
        public let rawValue: Int
        public init(rawValue: Int) { self.rawValue = rawValue }
        public static let x = Directions(rawValue: 1)
        public static let y = Directions(rawValue: 2)
        public static let z = Directions(rawValue: 4)
        public static let all: Directions = [.x, .y, .z]
    }

    public let structure: AtomicStructure
    /// The displacement of every atom [Å], one array per configuration (in the structure's atom order).
    public let displacements: [[SIMD3<Double>]]

    /// Given displacements (e.g. abTEM's, or from another generator).
    public init(structure: AtomicStructure, displacements: [[SIMD3<Double>]]) throws {
        guard !displacements.isEmpty else { throw SimulatoError.invalidFrozenPhonons("no configuration") }
        for (c, d) in displacements.enumerated() {
            guard d.count == structure.atoms.count else {
                throw SimulatoError.invalidFrozenPhonons("configuration \(c) has \(d.count) displacements for \(structure.atoms.count) atoms")
            }
            guard d.allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }) else {
                throw SimulatoError.invalidFrozenPhonons("configuration \(c) has a non-finite displacement")
            }
        }
        self.structure = structure
        self.displacements = displacements
    }

    /// Displacements drawn by the kit's generator: for configuration c the generator seeded with
    /// `configurationSeed(seed, c)`, three
    /// normals per atom in atom order (x, y, z), each axis not in `directions` left at zero, scaled by σ of the atom's Z
    /// (Å; every element of the structure must have a σ, as abTEM requires; 0 leaves it still).
    public init(structure: AtomicStructure, sigmas: [Int: Double], configurations: Int, seed: UInt64,
                directions: Directions = .all) throws {
        guard configurations >= 1 else { throw SimulatoError.invalidFrozenPhonons("configurations must be at least 1") }
        for z in Set(structure.atoms.map(\.atomicNumber)) where sigmas[z] == nil {
            throw SimulatoError.invalidFrozenPhonons("σ must be given for every element; Z = \(z) has none")
        }
        for (z, s) in sigmas where !(s.isFinite && s >= 0) {
            throw SimulatoError.invalidFrozenPhonons("σ of Z = \(z) must be finite and non-negative, got \(s)")
        }
        var all: [[SIMD3<Double>]] = []
        for c in 0..<configurations {
            var g = GaussianGenerator(seed: Self.configurationSeed(seed, c))
            all.append(structure.atoms.map { atom in
                let r = SIMD3(g.nextNormal(), g.nextNormal(), g.nextNormal())
                let s = sigmas[atom.atomicNumber] ?? 0
                return SIMD3(directions.contains(.x) ? s * r.x : 0, directions.contains(.y) ? s * r.y : 0,
                             directions.contains(.z) ? s * r.z : 0)
            })
        }
        try self.init(structure: structure, displacements: all)
    }

    public var count: Int { displacements.count }

    /// The generator seed of configuration c: the SplitMix64 finaliser of the finalised seed plus c, so the
    /// configurations of consecutive seeds do not coincide (seed + c made configuration c of seed s equal configuration
    /// c − 1 of seed s + 1).
    public static func configurationSeed(_ seed: UInt64, _ c: Int) -> UInt64 {
        func mix(_ v: UInt64) -> UInt64 {
            var z = v &+ 0x9E37_79B9_7F4A_7C15
            z = (z ^ (z >> 30)) &* 0xBF58_476D_1CE4_E5B9
            z = (z ^ (z >> 27)) &* 0x94D0_49BB_1331_11EB
            return z ^ (z >> 31)
        }
        return mix(mix(seed) &+ UInt64(c))
    }

    /// Configuration c: every atom at its position plus its displacement (not wrapped; the potential wraps it).
    public func configuration(_ c: Int) throws -> AtomicStructure {
        try AtomicStructure(box: structure.box, atoms: zip(structure.atoms, displacements[c]).map {
            AtomicStructure.Atom(atomicNumber: $0.0.atomicNumber, position: $0.0.position + $0.1)
        })
    }

    /// The projected potential of every configuration.
    public func potentials(nx: Int, ny: Int, sliceThickness: Double,
                           parameterization: LobatoParameterization = LobatoParameterization()) throws -> [ProjectedPotential] {
        try (0..<count).map {
            try ProjectedPotential(structure: configuration($0), nx: nx, ny: ny, sliceThickness: sliceThickness,
                                   parameterization: parameterization)
        }
    }

    /// The potential of every configuration, unbuilt (for `.metal`, which builds each slice on the GPU).
    public func lazyPotentials(nx: Int, ny: Int, sliceThickness: Double,
                               parameterization: LobatoParameterization = LobatoParameterization()) throws -> [Potential] {
        try (0..<count).map {
            try Potential(structure: configuration($0), nx: nx, ny: ny, sliceThickness: sliceThickness,
                          parameterization: parameterization)
        }
    }
}
