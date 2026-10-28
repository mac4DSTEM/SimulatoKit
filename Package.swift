// swift-tools-version: 6.2
// SimulatoKit — an electron-scattering simulation library: projected potentials, multislice propagation and frozen
// phonons (SimulatoMultislice), probes, aberrations and scans (SimulatoProbe), detectors and the STEM driver
// (SimulatoDetector) on a shared core of beams, grids, structures and waves (SimulatoCore). A Swift port of abTEM
// (GPL-3.0-or-later), checked against abTEM 1.0.10 (164e644) in float64 (docs/parity.md).
//
// Every accelerated path keeps a CPU reference path with a tested tolerance: the CPU path is float64 on Accelerate
// (vDSP), the Metal path float32 on MPSGraph, one API with the backend as a parameter.
//
// Swift 6 language mode, nonisolated by default, three upcoming features, no unsafe flags (a versioned dependency
// cannot carry them). Run the tests with `swift test -c release`. One product per module.
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
        // Constants, the electron beam, structures (a minimal local type, docs/adr/0001), grids, complex
        // fields, waves with their intensity and diffraction pattern, the antialias aperture, the float64 CPU FFT.
        .target(name: "SimulatoCore", swiftSettings: settings),
        // The Lobato parameterisation, slicing, the projected potential (infinite projection), transmission,
        // Fresnel propagation and the multislice loop on two backends (CPU float64, Metal float32).
        .target(name: "SimulatoMultislice", dependencies: ["SimulatoCore"], swiftSettings: settings),
        // The soft aperture, the polar aberrations C10 ... C56 and the probe at a position, grid and line scans, the CTF
        // with its envelopes and the closed forms, applied to waves for HRTEM images.
        .target(name: "SimulatoProbe", dependencies: ["SimulatoCore"], swiftSettings: settings),
        // Annular (BF, ADF, HAADF), pixelated, segmented, centre-of-mass and waves detectors and the STEM driver, averaged
        // over frozen-phonon configurations; line profiles and integrated centre of mass.
        .target(name: "SimulatoDetector", dependencies: ["SimulatoCore", "SimulatoMultislice", "SimulatoProbe"],
                swiftSettings: settings),
        // The first ten minutes: a CIF in, an HRTEM image and a STEM HAADF scan out (not a product).
        .executableTarget(name: "SimulatoExample",
                          dependencies: ["SimulatoCore", "SimulatoMultislice", "SimulatoProbe", "SimulatoDetector"],
                          path: "Examples/SimulatoExample", exclude: ["Al.cif"], swiftSettings: settings),
        // Test support both test sides use (relative-difference helpers); not a product.
        .target(name: "SimulatoKitTestSupport", dependencies: ["SimulatoCore"], path: "Tests/SimulatoKitTestSupport",
                swiftSettings: settings),
        .testTarget(name: "SimulatoCoreTests", dependencies: ["SimulatoCore"], swiftSettings: settings),
        .testTarget(name: "SimulatoMultisliceTests",
                    dependencies: ["SimulatoMultislice", "SimulatoProbe", "SimulatoCore", "SimulatoKitTestSupport"],
                    swiftSettings: settings),
        .testTarget(name: "SimulatoProbeTests", dependencies: ["SimulatoProbe", "SimulatoMultislice", "SimulatoCore"],
                    swiftSettings: settings),
        // The pinned comparisons against abTEM 1.0.10 (fixtures, truth, Metal against CPU, benchmarks): the slow half,
        // skipped with `--skip SimulatoKitParityTests` for a quick run, run by the plain `swift test`.
        .testTarget(name: "SimulatoKitParityTests",
                    dependencies: ["SimulatoCore", "SimulatoMultislice", "SimulatoProbe", "SimulatoDetector",
                                   "SimulatoKitTestSupport"],
                    exclude: ["Fixtures"], swiftSettings: settings),
        .testTarget(name: "SimulatoExampleTests", dependencies: ["SimulatoExample", "SimulatoCore"], swiftSettings: settings),
    ],
    swiftLanguageModes: [.v6]
)
