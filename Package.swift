// swift-tools-version: 6.2
// SimulatoKit — an electron-scattering simulation library: projected potentials and multislice
// propagation and frozen phonons (SimulatoMultislice), probes, aberrations and scans (SimulatoProbe), detectors and the
// STEM driver (SimulatoDetector) on a
// shared core of beams, grids, structures and waves (SimulatoCore). A Swift port of abTEM (GPL-3.0-or-later), checked
// against abTEM 1.0.10 (164e644) in float64; docs/registration-slice-1.md is the first registration.
//
// Every accelerated path keeps a CPU reference path with a registered tolerance: the CPU path is float64 on Accelerate
// (vDSP), the Metal path float32 on MPSGraph, one API with the backend as a parameter.
//
// Swift 6 language mode, nonisolated by default, three upcoming features, no unsafe flags
// (a versioned dependency cannot carry them): tools/test.sh runs the tests with `-c release`. One product per module.
import PackageDescription

let settings: [SwiftSetting] = [
    .enableUpcomingFeature("NonisolatedNonsendingByDefault"),
    .enableUpcomingFeature("InferIsolatedConformances"),
    .enableUpcomingFeature("MemberImportVisibility"),
]

let package = Package(
    name: "SimulatoKit",
    platforms: [.macOS(.v15), .iOS(.v18)],
    products: [
        .library(name: "SimulatoCore", targets: ["SimulatoCore"]),
        .library(name: "SimulatoMultislice", targets: ["SimulatoMultislice"]),
        .library(name: "SimulatoProbe", targets: ["SimulatoProbe"]),
        .library(name: "SimulatoDetector", targets: ["SimulatoDetector"]),
    ],
    targets: [
        // Constants, the electron beam, structures (a minimal local type until DiffractoKit's first tag), grids, complex
        // fields, waves with their intensity and diffraction pattern, the antialias aperture, the float64 CPU FFT.
        .target(name: "SimulatoCore", swiftSettings: settings),
        // Slice 1: the Lobato parameterisation, slicing, the projected potential (infinite projection), transmission,
        // Fresnel propagation and the multislice loop on two backends (CPU float64, Metal float32).
        .target(name: "SimulatoMultislice", dependencies: ["SimulatoCore"], swiftSettings: settings),
        // Slice 1: the soft aperture, the polar aberrations C10 ... C56 and the probe at a position; slice 2: scan
        // grids; slice 3: the CTF with its envelopes and the closed forms, applied to waves for HRTEM images.
        .target(name: "SimulatoProbe", dependencies: ["SimulatoCore"], swiftSettings: settings),
        // Slice 2: annular (BF, ADF, HAADF), pixelated and segmented detectors and the STEM driver, averaged over
        // frozen-phonon configurations (docs/registration-slice-2.md).
        .target(name: "SimulatoDetector", dependencies: ["SimulatoCore", "SimulatoMultislice", "SimulatoProbe"],
                swiftSettings: settings),
        // The first ten minutes: a CIF in, an HRTEM image and a STEM HAADF scan out (not a product).
        .executableTarget(name: "SimulatoExample",
                          dependencies: ["SimulatoCore", "SimulatoMultislice", "SimulatoProbe", "SimulatoDetector"],
                          path: "Examples/SimulatoExample", exclude: ["Al.cif"], swiftSettings: settings),
        .testTarget(name: "SimulatoCoreTests", dependencies: ["SimulatoCore"], swiftSettings: settings),
        .testTarget(name: "SimulatoMultisliceTests", dependencies: ["SimulatoMultislice", "SimulatoProbe"],
                    exclude: ["Fixtures"], swiftSettings: settings),
        .testTarget(name: "SimulatoProbeTests", dependencies: ["SimulatoProbe", "SimulatoMultislice", "SimulatoCore"],
                    exclude: ["Fixtures"], swiftSettings: settings),
        .testTarget(name: "SimulatoDetectorTests",
                    dependencies: ["SimulatoDetector", "SimulatoMultislice", "SimulatoProbe", "SimulatoCore"],
                    exclude: ["Fixtures"], swiftSettings: settings),
        .testTarget(name: "SimulatoExampleTests", dependencies: ["SimulatoExample", "SimulatoCore"], swiftSettings: settings),
    ],
    swiftLanguageModes: [.v6]
)
