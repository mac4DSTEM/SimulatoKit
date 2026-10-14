//
//  PhysicalConstants.swift
//  Role: the CODATA 2014 constants abTEM 1.0.10 reads from ASE 3.29.0 (ase.units), and ASE's unit factors, so the
//        wavelength and the interaction parameter equal abTEM's to the last bit (registration section 3.1).
//

/// CODATA 2014 values (SI), as ASE 3.29.0 and therefore abTEM 1.0.10 use them.
public enum PhysicalConstants {
    /// Planck constant h [J s].
    public static let planck = 6.62607004e-34
    /// Electron rest mass mₑ [kg].
    public static let electronMass = 9.10938356e-31
    /// Speed of light c [m/s].
    public static let speedOfLight = 299_792_458.0
    /// Elementary charge e [C].
    public static let elementaryCharge = 1.6021766208e-19
    /// Vacuum permittivity ε₀ [F/m].
    public static let vacuumPermittivity = 8.85418781762039e-12
    /// Bohr radius a₀ [Å].
    public static let bohrRadius = 0.5291772105638411
    /// abTEM's conversion factor κ = 4π ε₀ / (2π a₀ e) in ASE units (`abtem/core/constants.py`), the double abTEM computes.
    public static let kappa = 0.0208865737082965
}

/// ASE 3.29.0's unit factors (ASE units per SI unit), used only to reproduce abTEM's evaluation order.
package enum ASEUnits {
    package static let kg = 6.0221408585491615e+26
    package static let coulomb = 6.241509125883258e+18
    package static let second = 98226947884640.62
    package static let joule = 6.241509125883258e+18
}
