// swift-tools-version: 6.2
// SimulatoKit consumer check: a SEPARATE package that depends on SimulatoKit by path and uses, from outside the package,
// the whole public API of SimulatoCore, SimulatoMultislice, SimulatoProbe and SimulatoDetector. `package` access does not
// cross a package boundary, so this compiles only while each symbol it names is `public`. See README.md.
import Foundation
import PackageDescription

// A path dependency's identity is its directory's name (lower-cased): "simulatokit" in the repository, but not in a
// scratch copy; `.product` needs that name, so it is taken from the directory two levels above this manifest.
let kit = URL(fileURLWithPath: Context.packageDirectory).deletingLastPathComponent().deletingLastPathComponent()
    .lastPathComponent

let package = Package(
    name: "ConsumerCheck",
    platforms: [.macOS(.v15)],
    dependencies: [
        .package(path: "../.."),
    ],
    targets: [
        .executableTarget(name: "consumer-check", dependencies: [
            .product(name: "SimulatoCore", package: kit),
            .product(name: "SimulatoMultislice", package: kit),
            .product(name: "SimulatoProbe", package: kit),
            .product(name: "SimulatoDetector", package: kit),
        ]),
    ],
    swiftLanguageModes: [.v6]
)
