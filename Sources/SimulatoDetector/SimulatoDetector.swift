//
//  SimulatoDetector.swift
//  Role: the module's overview for DocC.
//

/// SimulatoDetector: the detectors of slice 2 (`AnnularDetector` for BF, ADF and HAADF, `PixelatedDetector` for 4D-STEM,
/// `SegmentedDetector` for DPC) and the STEM driver (`STEM.scan`), averaged over frozen-phonon configurations
/// (docs/registration-slice-2.md).
public enum SimulatoDetector {
    /// The scope of this module (slice 2).
    public static let plannedScope =
        "bright-field, annular (ADF, HAADF), pixelated (4D) and segmented detectors with frozen phonons (slice 2)"
}
