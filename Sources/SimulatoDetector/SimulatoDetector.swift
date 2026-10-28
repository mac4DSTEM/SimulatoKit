//
//  SimulatoDetector.swift
//  Role: the module's overview for DocC.
//

/// SimulatoDetector: the detectors (`AnnularDetector` for BF, ADF and HAADF, `PixelatedDetector` for 4D-STEM,
/// `SegmentedDetector` for DPC) and the STEM driver (`STEM.scan`), averaged over frozen-phonon configurations.
public enum SimulatoDetector {
    /// The scope of this module.
    public static let plannedScope =
        "bright-field, annular (ADF, HAADF), pixelated (4D) and segmented detectors with frozen phonons"
}
