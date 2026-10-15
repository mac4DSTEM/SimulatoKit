//
//  main.swift
//  Role: SimulatoKit's first ten minutes: a CIF in, two images out, through the public API only.
//        (a) HRTEM: a plane wave through the potential (Multislice.run), then the contrast transfer function
//            at Scherzer defocus (300 keV, Cs = 1 mm), written as <name>-hrtem.png;
//        (b) STEM HAADF: a probe scanned over one unit cell (STEM.scan, annular detector 70 to 200 mrad), written as
//            <name>-haadf.png; both on Metal when it is available (the potential built on the GPU), else on the CPU;
//        (c) every number printed and written to <name>-numbers.csv.
//
//          swift run -c release SimulatoExample [file.cif] [output-directory]
//
//        Without arguments: Examples/SimulatoExample/Al.cif (fcc Al, a = 4.0495 Å), written to the current directory.
//        The PNGs are 8-bit greyscale, stretched from the image's minimum to its maximum.
//

import CoreGraphics
import Foundation
import ImageIO
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe
import UniformTypeIdentifiers

/// An 8-bit greyscale PNG of `values` (x-major: index ix * ny + iy), stretched from min to max; x to the right, y down.
func writePNG(_ values: [Double], nx: Int, ny: Int, to url: URL) throws {
    let lo = values.min() ?? 0, hi = values.max() ?? 1
    var bytes = [UInt8](repeating: 0, count: nx * ny)
    for ix in 0..<nx {
        for iy in 0..<ny {
            let t = hi > lo ? (values[ix * ny + iy] - lo) / (hi - lo) : 0
            bytes[iy * nx + ix] = UInt8(max(0, min(255, (t * 255).rounded())))
        }
    }
    guard let provider = CGDataProvider(data: Data(bytes) as CFData),
          let image = CGImage(width: nx, height: ny, bitsPerComponent: 8, bitsPerPixel: 8, bytesPerRow: nx,
                              space: CGColorSpaceCreateDeviceGray(), bitmapInfo: CGBitmapInfo(rawValue: 0),
                              provider: provider, decode: nil, shouldInterpolate: false, intent: .defaultIntent),
          let dest = CGImageDestinationCreateWithURL(url as CFURL, UTType.png.identifier as CFString, 1, nil)
    else { throw MinimalCIF.Failure(description: "cannot write \(url.path)") }
    CGImageDestinationAddImage(dest, image, nil)
    guard CGImageDestinationFinalize(dest) else { throw MinimalCIF.Failure(description: "cannot write \(url.path)") }
}

func seconds(since t: Date) -> Double { Date().timeIntervalSince(t) }

do {
    // Arguments: an optional CIF path (ends in .cif) and an optional output directory.
    var cif = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Al.cif")
    var out = URL(fileURLWithPath: FileManager.default.currentDirectoryPath)
    for a in CommandLine.arguments.dropFirst() {
        if a.lowercased().hasSuffix(".cif") { cif = URL(fileURLWithPath: a) } else { out = URL(fileURLWithPath: a) }
    }
    try FileManager.default.createDirectory(at: out, withIntermediateDirectories: true)
    let (name, cell) = try MinimalCIF.read(String(contentsOf: cif, encoding: .utf8))

    // The crystal: 3 x 3 cells across, 10 cells deep, on a 256 x 256 grid, slices of about 2 Å, at 300 keV.
    let energy = 300e3, Cs = 1e7  // Å (1 mm)
    let structure = cell.repeated(3, 3, 10)
    let beam = try ElectronBeam(energy: energy)
    // A `Potential` is the sliced structure; on Metal each slice is built on the GPU, on the CPU by `build()`.
    let backend: Backend = Multislice.metalIsAvailable ? .metal : .cpu
    let t0 = Date()
    let potential = try Potential(structure: structure, nx: 256, ny: 256, sliceThickness: 2.0)
    let tPotential = seconds(since: t0)
    var rows: [(String, String)] = [
        ("structure", name), ("atoms", "\(structure.atoms.count)"), ("energy_keV", "\(energy / 1e3)"),
        ("wavelength_angstrom", "\(beam.wavelength)"), ("box_angstrom", "\(structure.box.x) \(structure.box.y) \(structure.box.z)"),
        ("gpts", "\(potential.grid.nx)x\(potential.grid.ny)"), ("slices", "\(potential.sliceCount)"), ("backend", backend == .metal ? "metal" : "cpu"),
        ("potential_seconds", "\(tPotential)"),
    ]

    // (a) HRTEM: plane wave -> multislice -> CTF at Scherzer defocus; the image is the intensity.
    let t1 = Date()
    let exit = try Multislice.run(.planeWave(beam: beam, grid: potential.grid), through: potential, backend: backend)[0]
    var aberrations = Aberrations()
    aberrations.Cs = Cs
    aberrations.defocus = try ContrastTransferFunction.scherzerDefocus(Cs: Cs, energy: energy)
    let ctf = try ContrastTransferFunction(aberrations: aberrations, focalSpread: 30)  // Δ = 30 Å
    let hrtem = try exit.applying(ctf).intensity()
    let mean = hrtem.sum / Double(hrtem.values.count)
    try writePNG(hrtem.values, nx: hrtem.nx, ny: hrtem.ny, to: out.appendingPathComponent("\(name)-hrtem.png"))
    rows += [("hrtem_seconds", "\(seconds(since: t1))"), ("scherzer_defocus_angstrom", "\(aberrations.defocus)"),
             ("point_resolution_angstrom", "\(try ContrastTransferFunction.pointResolution(Cs: Cs, energy: energy))"),
             ("hrtem_min", "\(hrtem.values.min()!)"), ("hrtem_max", "\(hrtem.max)"), ("hrtem_mean", "\(mean)")]
    print(String(format: "HRTEM  min %.6f  max %.6f  mean %.6f  (defocus %.1f Å)", hrtem.values.min()!, hrtem.max, mean,
                 aberrations.defocus))

    // (b) STEM HAADF: a 20 mrad probe over one unit cell (16 x 16 positions), detector 70 to 200 mrad.
    let t2 = Date()
    let probe = try Probe(beam: beam, grid: potential.grid, aperture: Aperture(semiangleCutoff: 20))
    let scan = try GridScan(start: SIMD2(0, 0), end: SIMD2(cell.box.x, cell.box.y), nx: 16, ny: 16)
    let haadf = try STEM.scan(probe: probe, potentials: [potential], scan: scan,
                              detectors: [AnnularDetector(inner: 70, outer: 200)], backend: backend)[0]
    try writePNG(haadf.values, nx: scan.nx, ny: scan.ny, to: out.appendingPathComponent("\(name)-haadf.png"))
    rows += [("stem_backend", backend == .metal ? "metal" : "cpu"), ("stem_seconds", "\(seconds(since: t2))"),
             ("haadf_scan", "\(scan.nx)x\(scan.ny)"), ("haadf_min", "\(haadf.values.min()!)"),
             ("haadf_max", "\(haadf.values.max()!)"), ("haadf_mean", "\(haadf.values.reduce(0, +) / Double(haadf.values.count))")]
    print("STEM HAADF (\(backend == .metal ? "metal" : "cpu"), rows are scan x, columns scan y), integrated intensity per probe:")
    for ix in 0..<scan.nx {
        print((0..<scan.ny).map { String(format: "%.5f", haadf.values[ix * scan.ny + $0]) }.joined(separator: " "))
    }

    // (c) the numbers, next to the PNGs.
    let csv = "quantity,value\n" + rows.map { "\($0.0),\($0.1)" }.joined(separator: "\n")
        + "\nhaadf_values_row_major_ix_iy,\(haadf.values.map { "\($0)" }.joined(separator: " "))\n"
    try csv.write(to: out.appendingPathComponent("\(name)-numbers.csv"), atomically: true, encoding: .utf8)
    print("total \(String(format: "%.2f", seconds(since: t0))) s; wrote \(name)-hrtem.png, \(name)-haadf.png, \(name)-numbers.csv in \(out.path)")
} catch {
    FileHandle.standardError.write(Data("SimulatoExample: \(error)\n".utf8))
    exit(1)
}
