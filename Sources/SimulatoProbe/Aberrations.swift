//
//  Aberrations.swift
//  Role: abTEM's polar aberration coefficients and the phase χ(α, φ) (registration section 3.8; `transfer.py`
//        `Aberrations._evaluate_from_angular_grid`): orders 1 to 5, coefficients in Å, angles in rad, defocus = −C10.
//

import Foundation
import SimulatoCore

/// Polar aberration coefficients (Krivanek's notation, abTEM's symbols). Lengths in Å, angles in rad.
public struct Aberrations: Sendable, Hashable {
    public var C10 = 0.0, C12 = 0.0, phi12 = 0.0
    public var C21 = 0.0, phi21 = 0.0, C23 = 0.0, phi23 = 0.0
    public var C30 = 0.0, C32 = 0.0, phi32 = 0.0, C34 = 0.0, phi34 = 0.0
    public var C41 = 0.0, phi41 = 0.0, C43 = 0.0, phi43 = 0.0, C45 = 0.0, phi45 = 0.0
    public var C50 = 0.0, C52 = 0.0, phi52 = 0.0, C54 = 0.0, phi54 = 0.0, C56 = 0.0, phi56 = 0.0

    public init() {}

    /// Defocus [Å], abTEM's sign: defocus = −C10.
    public var defocus: Double {
        get { -C10 }
        set { C10 = -newValue }
    }

    /// Spherical aberration Cs = C30 [Å].
    public var Cs: Double {
        get { C30 }
        set { C30 = newValue }
    }

    /// Every coefficient finite.
    var isFinite: Bool {
        [C10, C12, phi12, C21, phi21, C23, phi23, C30, C32, phi32, C34, phi34, C41, phi41, C43, phi43, C45, phi45,
         C50, C52, phi52, C54, phi54, C56, phi56].allSatisfy(\.isFinite)
    }

    /// χ [rad] at polar angle α [rad] and azimuth φ [rad] for wavelength λ [Å]:
    /// (2π/λ) Σ_n α^{n+1}/(n+1) Σ_m C_nm cos(m(φ − φ_nm)).
    public func chi(alpha: Double, phi: Double, wavelength: Double) -> Double {
        var x = 0.0
        if C10 != 0 || C12 != 0 || phi12 != 0 {
            x = x + 1 / 2 * pow(alpha, 2) * (C10 + C12 * cos(2 * (phi - phi12)))
        }
        if C21 != 0 || phi21 != 0 || C23 != 0 || phi23 != 0 {
            x = x + 1 / 3 * pow(alpha, 3) * (C21 * cos(phi - phi21) + C23 * cos(3 * (phi - phi23)))
        }
        if C30 != 0 || C32 != 0 || phi32 != 0 || C34 != 0 || phi34 != 0 {
            x = x + 1 / 4 * pow(alpha, 4) * (C30 + C32 * cos(2 * (phi - phi32)) + C34 * cos(4 * (phi - phi34)))
        }
        if C41 != 0 || phi41 != 0 || C43 != 0 || phi43 != 0 || C45 != 0 || phi45 != 0 {
            x = x + 1 / 5 * pow(alpha, 5)
                * (C41 * cos(phi - phi41) + C43 * cos(3 * (phi - phi43)) + C45 * cos(5 * (phi - phi45)))
        }
        if C50 != 0 || C52 != 0 || phi52 != 0 || C54 != 0 || phi54 != 0 || C56 != 0 || phi56 != 0 {
            x = x + 1 / 6 * pow(alpha, 6) * (C50 + C52 * cos(2 * (phi - phi52)) + C54 * cos(4 * (phi - phi54))
                                              + C56 * cos(6 * (phi - phi56)))
        }
        return x * (2 * Double.pi / wavelength)
    }
}
