//
//  SimulatoDetector.swift
//  Role: the module's placeholder until slice 2 is registered. Planned scope (suite ROADMAP.md, "Kits"; abTEM
//        `detectors.py`): the bright-field and annular detectors (ADF, HAADF; inner and outer angles in mrad), the
//        pixelated detector (a 4D-STEM cube per scan), the segmented detector (DPC), integrated over a scan and over
//        frozen-phonon configurations; checked against abTEM at the same pin. Depends on SimulatoCore.
//

/// The module's scope statement; replaced by the slice-2 API when it is registered.
public enum SimulatoDetector {
    /// What slice 2 will add here (docs/registration-slice-2.md, not yet written).
    public static let plannedScope =
        "bright-field, annular (ADF, HAADF), pixelated (4D) and segmented detectors with frozen phonons (slice 2)"
}
