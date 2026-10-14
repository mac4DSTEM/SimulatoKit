//
//  main.swift
//  Role: SimulatoKit's first ten minutes: read a CIF (orthogonal cell), repeat it, build the projected potential, run a
//        plane wave through it, and write the exit-wave intensity and the diffraction pattern as PNGs plus the numbers as
//        JSON. Usage:
//
//          swift run -c release SimulatoExample [file.cif] [--energy keV] [--repeat nx ny nz] [--gpts n]
//                                               [--slice Å] [--backend cpu|metal] [--out directory]
//
//        Without a file, it uses Examples/SimulatoExample/Al.cif (fcc Al, a = 4.0495 Å).
//

import CoreGraphics
import Foundation
import ImageIO
import SimulatoCore
import SimulatoMultislice
import UniformTypeIdentifiers

struct Options {
    var cif = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Al.cif")
    var energy = 200.0
    var repeats = (4, 4, 10)
    var gpts = 256
    var slice = 2.0
    var backend = Backend.cpu
    var out = URL(fileURLWithPath: "Examples/output")
}

func parse() throws -> Options {
    var o = Options()
    var args = Array(CommandLine.arguments.dropFirst())
    func next() throws -> String {
        guard !args.isEmpty else { throw MinimalCIF.Failure(description: "an option is missing its value") }
        return args.removeFirst()
    }
    while !args.isEmpty {
        let a = args.removeFirst()
        switch a {
        case "--energy": o.energy = Double(try next()) ?? o.energy
        case "--repeat": o.repeats = (Int(try next()) ?? 1, Int(try next()) ?? 1, Int(try next()) ?? 1)
        case "--gpts": o.gpts = Int(try next()) ?? o.gpts
        case "--slice": o.slice = Double(try next()) ?? o.slice
        case "--backend": o.backend = try next() == "metal" ? .metal : .cpu
        case "--out": o.out = URL(fileURLWithPath: try next())
        default: o.cif = URL(fileURLWithPath: a)
        }
    }
    return o
}

/// An 8-bit greyscale PNG; pixel (ix, iy) of the field is column ix, row iy (x to the right, y down).
func writePNG(_ field: RealField, to url: URL, logScale: Bool) throws {
    var v = field.values
    if logScale { let floor = (v.max() ?? 1) * 1e-6; v = v.map { log10(max($0, floor)) } }
    let lo = v.min() ?? 0, hi = v.max() ?? 1
    var bytes = [UInt8](repeating: 0, count: field.nx * field.ny)
    for ix in 0..<field.nx {
        for iy in 0..<field.ny {
            let t = hi > lo ? (v[ix * field.ny + iy] - lo) / (hi - lo) : 0
            bytes[iy * field.nx + ix] = UInt8(max(0, min(255, (t * 255).rounded())))
        }
    }
    let provider = CGDataProvider(data: Data(bytes) as CFData)!
    let image = CGImage(width: field.nx, height: field.ny, bitsPerComponent: 8, bitsPerPixel: 8, bytesPerRow: field.nx,
                        space: CGColorSpaceCreateDeviceGray(), bitmapInfo: CGBitmapInfo(rawValue: 0), provider: provider,
                        decode: nil, shouldInterpolate: false, intent: .defaultIntent)!
    guard let dest = CGImageDestinationCreateWithURL(url as CFURL, UTType.png.identifier as CFString, 1, nil) else {
        throw MinimalCIF.Failure(description: "cannot write \(url.path)")
    }
    CGImageDestinationAddImage(dest, image, nil)
    guard CGImageDestinationFinalize(dest) else { throw MinimalCIF.Failure(description: "cannot write \(url.path)") }
}

do {
    let o = try parse()
    let (name, cell) = try MinimalCIF.read(String(contentsOf: o.cif, encoding: .utf8))
    let structure = cell.repeated(o.repeats.0, o.repeats.1, o.repeats.2)
    let beam = try ElectronBeam(energy: o.energy * 1e3)
    let start = Date()
    let potential = try ProjectedPotential(structure: structure, nx: o.gpts, ny: o.gpts, sliceThickness: o.slice)
    let exit = try Multislice.run(.planeWave(beam: beam, grid: potential.grid), through: potential, backend: o.backend)[0]
    let seconds = Date().timeIntervalSince(start)
    let image = exit.intensity()
    let pattern = try exit.diffractionPattern()
    try FileManager.default.createDirectory(at: o.out, withIntermediateDirectories: true)
    try writePNG(image, to: o.out.appendingPathComponent("\(name)-exit-intensity.png"), logScale: false)
    try writePNG(pattern, to: o.out.appendingPathComponent("\(name)-diffraction-log.png"), logScale: true)
    // The ten strongest beams, (kx, ky) in 1/Å from the pattern's centre.
    let cx = pattern.nx / 2, cy = pattern.ny / 2
    let beams = pattern.values.indices.sorted { pattern.values[$0] > pattern.values[$1] }.prefix(10).map { i -> [String: Double] in
        let ix = i / pattern.ny, iy = i % pattern.ny
        return ["kx": Double(ix - cx) * pattern.samplingX, "ky": Double(iy - cy) * pattern.samplingY,
                "intensity": pattern.values[i]]
    }
    let numbers: [String: Any] = [
        "simulatokit_reference": SimulatoCore.referenceVersion,
        "structure": name, "atoms": structure.atoms.count,
        "box_angstrom": [structure.box.x, structure.box.y, structure.box.z],
        "energy_eV": beam.energy, "wavelength_angstrom": beam.wavelength, "sigma_per_V_angstrom": beam.interactionParameter,
        "gpts": [potential.grid.nx, potential.grid.ny], "sampling_angstrom": [potential.grid.samplingX, potential.grid.samplingY],
        "slices": potential.slices.count, "slice_thickness_angstrom": potential.slices[0].thickness,
        "backend": o.backend == .metal ? "metal" : "cpu", "seconds": seconds,
        "exit_total_intensity": exit.totalIntensity, "exit_intensity_min": image.values.min()!, "exit_intensity_max": image.max,
        "pattern_gpts": [pattern.nx, pattern.ny], "pattern_sum": pattern.sum, "strongest_beams": beams,
    ]
    let json = try JSONSerialization.data(withJSONObject: numbers, options: [.prettyPrinted, .sortedKeys])
    try json.write(to: o.out.appendingPathComponent("\(name)-numbers.json"))
    print(String(data: json, encoding: .utf8)!)
    print("wrote \(o.out.path)/\(name)-exit-intensity.png, -diffraction-log.png, -numbers.json")
} catch {
    FileHandle.standardError.write(Data("SimulatoExample: \(error)\n".utf8))
    exit(1)
}
