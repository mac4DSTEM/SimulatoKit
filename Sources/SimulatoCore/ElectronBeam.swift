//
//  ElectronBeam.swift
//  Role: the beam energy and the two numbers multislice needs from it: the relativistic wavelength λ and the
//        interaction parameter σ (abTEM `core/energy.py`, same evaluation order).
//

import Foundation

/// An electron beam of a given kinetic energy.
public struct ElectronBeam: Sendable, Hashable {
    /// Kinetic energy [eV].
    public let energy: Double

    /// - Throws: `SimulatoError.invalidEnergy` unless the energy is finite and positive.
    public init(energy: Double) throws {
        guard energy.isFinite, energy > 0 else { throw SimulatoError.invalidEnergy(energy) }
        self.energy = energy
    }

    /// Relativistic de Broglie wavelength λ [Å]: h c / √(E (2 mₑ c² / e + E)) / e × 10¹⁰.
    public var wavelength: Double { Self.wavelength(energy: energy) }

    /// Interaction parameter σ [1/(V·Å)]: 2π m(E) e λ / h², m(E) = mₑ (1 + e E / (mₑ c²)).
    public var interactionParameter: Double {
        let c = PhysicalConstants.self
        let mass = (1 + c.elementaryCharge * energy / (c.electronMass * (c.speedOfLight * c.speedOfLight)))
            * c.electronMass
        let hsj = c.planck * ASEUnits.second * ASEUnits.joule
        return 2 * Double.pi * mass * ASEUnits.kg * c.elementaryCharge * ASEUnits.coulomb * wavelength / (hsj * hsj)
    }

    /// Relativistic mass correction factor 1 + e E / (mₑ c²).
    public var relativisticMassFactor: Double {
        let c = PhysicalConstants.self
        return 1 + c.elementaryCharge * energy / (c.electronMass * (c.speedOfLight * c.speedOfLight))
    }

    static func wavelength(energy: Double) -> Double {
        let c = PhysicalConstants.self
        let c2 = c.speedOfLight * c.speedOfLight
        return c.planck * c.speedOfLight / (energy * (2 * c.electronMass * c2 / c.elementaryCharge + energy)).squareRoot()
            / c.elementaryCharge * 1.0e10
    }
}
