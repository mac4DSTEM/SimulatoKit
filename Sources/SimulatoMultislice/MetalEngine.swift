//
//  MetalEngine.swift
//  Role: the Metal path on the GPU end to end (registration-performance.md section 1): one compiled MPSGraph executable
//        per step and shape, encoded once per slice into one command buffer and committed once (or at each requested
//        exit plane); the projected potential built in the step from host-merged bilinear deltas and per-element kernels
//        (for a `Potential`), or the transmission built from an uploaded float32 V (for a `ProjectedPotential`); the
//        wave kept on the device in two ping-pong buffers; probes and diffraction intensities for scans on the GPU.
//        Float32 throughout on the device; every array the step multiplies by (kernels, aperture, propagators, the
//        probe's spectrum and shift ramps) is computed on the host in float64 by the CPU path's code and cast once.
//

import Foundation
import Metal
import MetalPerformanceShaders
import MetalPerformanceShadersGraph
import Synchronization
import SimulatoCore

/// What the Metal path computes per slice: from atoms (a `Potential`) or from a built V (a `ProjectedPotential`).
enum MetalSliceSource: Hashable {
    case atoms(elements: Int, contributions: Int)
    case values
}

/// Options of the Metal path that tests compare (not public).
struct MetalOptions: Hashable {
    /// Rescale the transmission and the wave after each FFT pair so their Σ|·|² follows Parseval as measured on the
    /// device (registration-performance.md 1.7, amendment 1).
    var renormalize: Bool

    static let standard = MetalOptions(renormalize: true)
}

/// A compiled executable with its feeds in the executable's order.
final class MetalExecutable: @unchecked Sendable {
    let executable: MPSGraphExecutable
    let feeds: [MPSGraphTensor]
    let output: [Int]

    init(_ graph: MPSGraph, device: MPSGraphDevice, feeds: [(MPSGraphTensor, [Int], MPSDataType)],
         targets: [MPSGraphTensor]) {
        var shaped: [MPSGraphTensor: MPSGraphShapedType] = [:]
        for (t, shape, type) in feeds {
            shaped[t] = MPSGraphShapedType(shape: shape.map { NSNumber(value: $0) }, dataType: type)
        }
        executable = graph.compile(with: device, feeds: shaped, targetTensors: targets, targetOperations: nil,
                                   compilationDescriptor: nil)
        self.feeds = executable.feedTensors ?? feeds.map(\.0)
        output = []
    }

    /// The inputs in the executable's order.
    func inputs(_ map: [MPSGraphTensor: MPSGraphTensorData]) -> [MPSGraphTensorData] { feeds.map { map[$0]! } }
}

/// The device, its queue and the compiled executables, shared by every Metal run in the process.
final class MetalEngine: @unchecked Sendable {
    let device: any MTLDevice
    let queue: any MTLCommandQueue
    let graphDevice: MPSGraphDevice
    /// One Metal run at a time (the executables and the queue are shared).
    let runLock = NSLock()
    private let cache = Mutex<[String: MetalStep]>([:])
    /// The largest `currentAllocatedSize` seen at a commit since the last reset (registration-performance.md 1.8).
    let peak = Mutex<Int>(0)

    func notePeak() {
        let now = device.currentAllocatedSize
        peak.withLock { $0 = max($0, now) }
    }

    static let shared: MetalEngine? = MetalEngine()

    private init?() {
        guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else { return nil }
        self.device = device
        self.queue = queue
        graphDevice = MPSGraphDevice(mtlDevice: device)
    }

    static func require() throws -> MetalEngine {
        guard let e = shared else {
            throw SimulatoError.metalUnavailable("no Metal device (MTLCreateSystemDefaultDevice() returned nil)")
        }
        return e
    }

    /// The executable for `key`, built by `make` the first time.
    func step(_ key: String, _ make: () -> MetalStep) -> MetalStep {
        if let s = cache.withLock({ $0[key] }) { return s }
        let s = make()
        cache.withLock { $0[key] = s }
        return s
    }

    /// A shared-storage tensor of `count` 32-bit or 64-bit elements, zeroed.
    func tensor(_ shape: [Int], _ type: MPSDataType) throws -> (MPSGraphTensorData, any MTLBuffer) {
        let bytes = shape.reduce(1, *) * (type == .complexFloat32 ? 8 : 4)
        guard let buffer = device.makeBuffer(length: max(bytes, 4), options: .storageModeShared) else {
            throw SimulatoError.metalFailure("the device made no buffer of \(bytes) bytes")
        }
        memset(buffer.contents(), 0, max(bytes, 4))
        return (MPSGraphTensorData(buffer, shape: shape.map { NSNumber(value: $0) }, dataType: type), buffer)
    }

    func floats(_ values: [Float], _ shape: [Int], _ type: MPSDataType = .float32) throws -> MPSGraphTensorData {
        let (t, b) = try tensor(shape, type)
        values.withUnsafeBytes { b.contents().copyMemory(from: $0.baseAddress!, byteCount: $0.count) }
        return t
    }

    func ints(_ values: [Int32], _ shape: [Int]) throws -> MPSGraphTensorData {
        let (t, b) = try tensor(shape, .int32)
        values.withUnsafeBytes { b.contents().copyMemory(from: $0.baseAddress!, byteCount: $0.count) }
        return t
    }

    /// A complex field as complex64 (interleaved), cast from float64.
    func complex(_ f: ComplexField, shape: [Int]? = nil) throws -> MPSGraphTensorData {
        var v = [Float](repeating: 0, count: 2 * f.real.count)
        for i in f.real.indices {
            v[2 * i] = Float(f.real[i])
            v[2 * i + 1] = Float(f.imaginary[i])
        }
        return try floats(v, shape ?? [f.nx, f.ny], .complexFloat32)
    }

    func check(_ cb: MPSCommandBuffer, _ what: String) throws {
        if cb.status == .error {
            throw SimulatoError.metalFailure("\(what): \(cb.error.map { "\($0)" } ?? "the command buffer failed")")
        }
    }
}

/// One compiled graph with its named placeholders.
struct MetalStep: @unchecked Sendable {
    let executable: MetalExecutable
    let placeholders: [String: MPSGraphTensor]

    func inputs(_ named: [String: MPSGraphTensorData]) -> [MPSGraphTensorData] {
        var map: [MPSGraphTensor: MPSGraphTensorData] = [:]
        for (name, data) in named { map[placeholders[name]!] = data }
        return executable.inputs(map)
    }
}

/// The graph pieces (registration-performance.md 1.3–1.6).
enum MetalGraphs {
    static func shape(_ a: [Int]) -> [NSNumber] { a.map { NSNumber(value: $0) } }

    /// The forward FFT, unscaled (NumPy's `fft2`).
    static var forward: MPSGraphFFTDescriptor {
        let d = MPSGraphFFTDescriptor()
        d.inverse = false
        d.scalingMode = .none
        return d
    }

    /// The inverse FFT scaled by 1/N (NumPy's `ifft2`).
    static var inverse: MPSGraphFFTDescriptor {
        let d = MPSGraphFFTDescriptor()
        d.inverse = true
        d.scalingMode = .size
        return d
    }

    static func energy(_ g: MPSGraph, _ x: MPSGraphTensor, axes: [NSNumber]) -> MPSGraphTensor {
        let re = g.realPartOfTensor(tensor: x, name: nil), im = g.imaginaryPartOfTensor(tensor: x, name: nil)
        return g.reductionSum(with: g.addition(g.multiplication(re, re, name: nil), g.multiplication(im, im, name: nil),
                                               name: nil), axes: axes, name: nil)
    }

    /// y · √((E_Y E_x) / (E_X E_y)): the round trip's Σ|·|² made to follow Parseval as the device measured it before and
    /// after the forward transform (1 where a product is 0).
    static func renormalized(_ g: MPSGraph, _ y: MPSGraphTensor, ex: MPSGraphTensor, eX: MPSGraphTensor, eY: MPSGraphTensor,
                             axes: [NSNumber]) -> MPSGraphTensor {
        let ey = energy(g, y, axes: axes)
        let den = g.multiplication(eX, ey, name: nil)
        let num = g.multiplication(eY, ex, name: nil)
        let zero = g.constant(0, dataType: .float32), one = g.constant(1, dataType: .float32)
        let s = g.select(predicate: g.greaterThan(den, zero, name: nil),
                         trueTensor: g.squareRoot(with: g.division(num, den, name: nil), name: nil),
                         falseTensor: one, name: nil)
        return g.multiplication(y, g.cast(s, to: .complexFloat32, name: "s"), name: nil)
    }

    /// V [nx, ny] from the merged deltas of E elements: Re IFFT(Σ_e FFT(δ_e) · K_e).
    static func potential(_ g: MPSGraph, idx: MPSGraphTensor, w: MPSGraphTensor, kernels: MPSGraphTensor,
                          elements: Int, nx: Int, ny: Int) -> MPSGraphTensor {
        let d = g.scatter(w, indices: idx, shape: shape([elements * nx * ny]), axis: 0, mode: .add, name: nil)
        let d3 = g.cast(g.reshape(d, shape: shape([elements, nx, ny]), name: nil), to: .complexFloat32, name: "d")
        let f = g.fastFourierTransform(d3, axes: [1, 2], descriptor: forward, name: nil)
        let fk = g.multiplication(f, g.cast(kernels, to: .complexFloat32, name: "k"), name: nil)
        let sum = g.reductionSum(with: fk, axis: 0, name: nil)
        let v = g.realPartOfTensor(tensor: g.fastFourierTransform(sum, axes: [1, 2], descriptor: inverse, name: nil),
                                   name: nil)
        return g.reshape(v, shape: shape([nx, ny]), name: nil)
    }

    /// t = IFFT(FFT(exp(iσV)) · A), [nx, ny].
    static func transmission(_ g: MPSGraph, v: MPSGraphTensor, sigma: MPSGraphTensor, aperture: MPSGraphTensor,
                             renormalize: Bool) -> MPSGraphTensor {
        let phase = g.multiplication(v, sigma, name: nil)
        let e = g.complexTensor(realTensor: g.cos(with: phase, name: nil), imaginaryTensor: g.sin(with: phase, name: nil),
                                name: nil)
        let spectrum = g.fastFourierTransform(e, axes: [0, 1], descriptor: forward, name: nil)
        let limited = g.multiplication(spectrum, g.cast(aperture, to: .complexFloat32, name: "a"), name: nil)
        let t = g.fastFourierTransform(limited, axes: [0, 1], descriptor: inverse, name: nil)
        guard renormalize else { return t }
        let axes: [NSNumber] = [0, 1]
        return renormalized(g, t, ex: energy(g, e, axes: axes), eX: energy(g, spectrum, axes: axes),
                            eY: energy(g, limited, axes: axes), axes: axes)
    }

    /// ψ ← IFFT(FFT(ψ · t) · P) over a batch [B, nx, ny].
    static func propagate(_ g: MPSGraph, psi: MPSGraphTensor, t: MPSGraphTensor, p: MPSGraphTensor,
                          renormalize: Bool) -> MPSGraphTensor {
        let x = g.multiplication(psi, t, name: nil)
        let spectrum = g.fastFourierTransform(x, axes: [1, 2], descriptor: forward, name: nil)
        let propagated = g.multiplication(spectrum, p, name: nil)
        let y = g.fastFourierTransform(propagated, axes: [1, 2], descriptor: inverse, name: nil)
        guard renormalize else { return y }
        let axes: [NSNumber] = [1, 2]
        return renormalized(g, y, ex: energy(g, x, axes: axes), eX: energy(g, spectrum, axes: axes),
                            eY: energy(g, propagated, axes: axes), axes: axes)
    }
}

/// A potential prepared on the device for one beam: per-slice data, the kernels, the aperture, the propagators.
struct MetalPreparedPotential: @unchecked Sendable {
    let grid: Grid
    let source: MetalSliceSource
    /// `.atoms`: (indices, weights) per slice; `.values`: V per slice in `values`.
    let indices: [MPSGraphTensorData]
    let weights: [MPSGraphTensorData]
    let values: [MPSGraphTensorData]
    let kernels: MPSGraphTensorData?
    let aperture: MPSGraphTensorData
    let propagators: [MPSGraphTensorData]
    let sigma: MPSGraphTensorData
    var sliceCount: Int { propagators.count }
}

extension MetalMultislice {
    /// The contributions of slice j as unique (flat index e · N + ix · ny + iy, weight) pairs: the bilinear weights of
    /// `SliceBuilder.superposeDeltas`, summed per pixel in float64 on the host (registration-performance.md 1.3).
    static func contributions(_ potential: Potential, elements: [Int]) -> [[(Int, Double)]] {
        let grid = potential.grid
        let n = grid.count, nx = grid.nx, ny = grid.ny
        let element = Dictionary(uniqueKeysWithValues: elements.enumerated().map { ($0.element, $0.offset) })
        var perSlice = [[Int: Double]](repeating: [:], count: potential.sliceCount)
        func mod(_ a: Int, _ m: Int) -> Int { let r = a % m; return r < 0 ? r + m : r }
        for (a, atom) in potential.atoms.enumerated() {
            let j = potential.sliceOfAtom[a]
            let e = element[atom.atomicNumber]!
            let px = atom.position.x / grid.samplingX, py = atom.position.y / grid.samplingY
            let r = Int(px.rounded(.down)), c = Int(py.rounded(.down))
            let x = px - Double(r), y = py - Double(c)
            let xy = x * y
            let w = [1 + xy - y - x, x - xy, y - xy, xy]
            for k in 0..<4 {
                let i = mod(r + (k == 1 || k == 3 ? 1 : 0), nx)
                let jj = mod(c + (k >= 2 ? 1 : 0), ny)
                perSlice[j][e * n + i * ny + jj, default: 0] += w[k]
            }
        }
        return perSlice.map { $0.sorted { $0.key < $1.key }.map { ($0.key, $0.value) } }
    }

    static func prepare(_ potential: Potential, beam: ElectronBeam, engine: MetalEngine) throws -> MetalPreparedPotential {
        let grid = potential.grid
        let elements = Set(potential.atoms.map(\.atomicNumber)).sorted()
        let e = max(elements.count, 1)
        let lists = contributions(potential, elements: elements)
        // Padded to a power of two (at least 16) so potentials of similar size share one executable; padding adds 0 at 0.
        var m = 16
        while m < (lists.map(\.count).max() ?? 0) { m *= 2 }
        var indices: [MPSGraphTensorData] = [], weights: [MPSGraphTensorData] = []
        for list in lists {
            var idx = [Int32](repeating: 0, count: m), w = [Float](repeating: 0, count: m)
            for (k, (i, v)) in list.enumerated() {
                idx[k] = Int32(i)
                w[k] = Float(v)
            }
            indices.append(try engine.ints(idx, [m]))
            weights.append(try engine.floats(w, [m]))
        }
        let builder = try SliceBuilder(grid: grid, parameterization: potential.parameterization)
        var k = [Float](repeating: 0, count: e * grid.count)
        for (q, z) in elements.enumerated() {
            let kernel = try builder.kernel(z)
            for i in 0..<grid.count { k[q * grid.count + i] = Float(kernel[i]) }
        }
        return MetalPreparedPotential(
            grid: grid, source: .atoms(elements: e, contributions: m), indices: indices, weights: weights, values: [],
            kernels: try engine.floats(k, [e, grid.nx, grid.ny]),
            aperture: try engine.floats(Antialias.aperture(on: grid).map(Float.init), [grid.nx, grid.ny]),
            propagators: try propagators(potential.sliceThicknesses, grid: grid, beam: beam, engine: engine),
            sigma: try engine.floats([Float(beam.interactionParameter)], [1]))
    }

    static func prepare(_ potential: ProjectedPotential, beam: ElectronBeam,
                        engine: MetalEngine) throws -> MetalPreparedPotential {
        let grid = potential.grid
        return MetalPreparedPotential(
            grid: grid, source: .values, indices: [], weights: [],
            values: try potential.slices.map { try engine.floats($0.values.map(Float.init), [grid.nx, grid.ny]) },
            kernels: nil,
            aperture: try engine.floats(Antialias.aperture(on: grid).map(Float.init), [grid.nx, grid.ny]),
            propagators: try propagators(potential.slices.map(\.thickness), grid: grid, beam: beam, engine: engine),
            sigma: try engine.floats([Float(beam.interactionParameter)], [1]))
    }

    /// One complex64 propagator per distinct thickness, computed by the CPU path's code in float64.
    static func propagators(_ thicknesses: [Double], grid: Grid, beam: ElectronBeam,
                            engine: MetalEngine) throws -> [MPSGraphTensorData] {
        var byThickness: [Double: MPSGraphTensorData] = [:]
        return try thicknesses.map { dz in
            if let p = byThickness[dz] { return p }
            let p = try engine.complex(Multislice.propagator(grid: grid, beam: beam, thickness: dz))
            byThickness[dz] = p
            return p
        }
    }

    /// The step executable for a prepared potential and a batch of b waves.
    static func stepExecutable(_ p: MetalPreparedPotential, batch b: Int, options: MetalOptions,
                               engine: MetalEngine) -> MetalStep {
        let nx = p.grid.nx, ny = p.grid.ny
        return engine.step("step \(p.source) \(nx)x\(ny) b\(b) \(options)") {
            let g = MPSGraph()
            let psi = g.placeholder(shape: MetalGraphs.shape([b, nx, ny]), dataType: .complexFloat32, name: "psi")
            let aperture = g.placeholder(shape: MetalGraphs.shape([nx, ny]), dataType: .float32, name: "aperture")
            let prop = g.placeholder(shape: MetalGraphs.shape([nx, ny]), dataType: .complexFloat32, name: "p")
            let sigma = g.placeholder(shape: [1], dataType: .float32, name: "sigma")
            var feeds: [(MPSGraphTensor, [Int], MPSDataType)] = [
                (psi, [b, nx, ny], .complexFloat32), (aperture, [nx, ny], .float32), (prop, [nx, ny], .complexFloat32),
                (sigma, [1], .float32),
            ]
            var named = ["psi": psi, "aperture": aperture, "p": prop, "sigma": sigma]
            let v: MPSGraphTensor
            switch p.source {
            case .atoms(let e, let m):
                let idx = g.placeholder(shape: [NSNumber(value: m)], dataType: .int32, name: "idx")
                let w = g.placeholder(shape: [NSNumber(value: m)], dataType: .float32, name: "w")
                let k = g.placeholder(shape: MetalGraphs.shape([e, nx, ny]), dataType: .float32, name: "kernels")
                feeds += [(idx, [m], .int32), (w, [m], .float32), (k, [e, nx, ny], .float32)]
                named["idx"] = idx
                named["w"] = w
                named["kernels"] = k
                v = MetalGraphs.potential(g, idx: idx, w: w, kernels: k, elements: e, nx: nx, ny: ny)
            case .values:
                let values = g.placeholder(shape: MetalGraphs.shape([nx, ny]), dataType: .float32, name: "v")
                feeds.append((values, [nx, ny], .float32))
                named["v"] = values
                v = values
            }
            let t = MetalGraphs.transmission(g, v: v, sigma: sigma, aperture: aperture, renormalize: options.renormalize)
            let out = MetalGraphs.propagate(g, psi: psi, t: t, p: prop, renormalize: options.renormalize)
            return MetalStep(executable: MetalExecutable(g, device: engine.graphDevice, feeds: feeds, targets: [out]),
                             placeholders: named)
        }
    }

    static func stepInputs(_ p: MetalPreparedPotential, slice j: Int, psi: MPSGraphTensorData) -> [String: MPSGraphTensorData] {
        var named = ["psi": psi, "aperture": p.aperture, "p": p.propagators[j], "sigma": p.sigma]
        switch p.source {
        case .atoms:
            named["idx"] = p.indices[j]
            named["w"] = p.weights[j]
            named["kernels"] = p.kernels!
        case .values:
            named["v"] = p.values[j]
        }
        return named
    }

    /// Propagates the batch in `psi` (shape [b, nx, ny], overwritten) through every slice: one command buffer, committed
    /// at each plane in `planes` (zero-based slice indices; the last slice always) where `read` gets the wave. Every
    /// slice and every batch runs in its own autorelease pool: without one, each batch's buffers lived to the end of the
    /// scan (the first benchmark: 73 GB at 1024², 4096 positions).
    static func propagate(_ p: MetalPreparedPotential, psi: MPSGraphTensorData, batch b: Int, planes: Set<Int>,
                          options: MetalOptions, engine: MetalEngine,
                          read: (Int, MPSGraphTensorData) throws -> Void) throws {
        let step = stepExecutable(p, batch: b, options: options, engine: engine)
        var current = psi
        var other = try engine.tensor([b, p.grid.nx, p.grid.ny], .complexFloat32).0
        var cb = MPSCommandBuffer(from: engine.queue)
        let last = p.sliceCount - 1
        for j in 0...last {
            try autoreleasepool {
                _ = step.executable.executable.encode(to: cb, inputs: step.inputs(stepInputs(p, slice: j, psi: current)),
                                                      results: [other], executionDescriptor: nil)
                swap(&current, &other)
                if planes.contains(j) || j == last {
                    engine.notePeak()
                    cb.commit()
                    cb.waitUntilCompleted()
                    try engine.check(cb, "slice \(j)")
                    try read(j, current)
                    if j < last { cb = MPSCommandBuffer(from: engine.queue) }
                }
            }
        }
    }

    /// Complex64 [b, nx, ny] as b float64 fields.
    static func fields(_ data: MPSGraphTensorData, batch b: Int, grid: Grid) -> [ComplexField] {
        var interleaved = [Float](repeating: 0, count: 2 * grid.count * b)
        data.mpsndarray().readBytes(&interleaved, strideBytes: nil)
        return (0..<b).map { k in
            var f = ComplexField(nx: grid.nx, ny: grid.ny)
            let base = 2 * grid.count * k
            for i in 0..<grid.count {
                f.real[i] = Double(interleaved[base + 2 * i])
                f.imaginary[i] = Double(interleaved[base + 2 * i + 1])
            }
            return f
        }
    }

    /// The batch size whose waves fit in 512 MB (registration-slice-2.md 3.6).
    static func batchSize(_ grid: Grid) -> Int { max(1, (512 << 20) / (grid.count * 8 * 3)) }

    /// The waves of one incident at the requested planes (−1: the incident itself, which the caller adds).
    static func run(_ incident: Waves, prepared p: MetalPreparedPotential, planes: [Int],
                    options: MetalOptions = .standard) throws -> [Waves] {
        let engine = try MetalEngine.require()
        engine.runLock.lock()
        defer { engine.runLock.unlock() }
        var result: [Waves] = []
        if planes.first == -1 { result.append(incident) }
        let wanted = Set(planes.filter { $0 >= 0 })
        let psi = try engine.complex(incident.field, shape: [1, p.grid.nx, p.grid.ny])
        try propagate(p, psi: psi, batch: 1, planes: wanted, options: options, engine: engine) { j, data in
            guard wanted.contains(j) else { return }
            result.append(Waves(field: fields(data, batch: 1, grid: p.grid)[0], grid: p.grid, beam: incident.beam,
                                normalization: incident.normalization))
        }
        return result
    }

    /// The exit waves of a batch of incidents (after the last slice), in input order, in batches that fit 512 MB.
    static func runBatch(_ incidents: [ComplexField], prepared p: MetalPreparedPotential,
                         options: MetalOptions = .standard) throws -> [ComplexField] {
        let engine = try MetalEngine.require()
        engine.runLock.lock()
        defer { engine.runLock.unlock() }
        let size = batchSize(p.grid)
        var out: [ComplexField] = []
        for start in stride(from: 0, to: incidents.count, by: size) {
            try autoreleasepool {
                let chunk = Array(incidents[start..<min(start + size, incidents.count)])
                var v = [Float](repeating: 0, count: 2 * p.grid.count * chunk.count)
                for (k, f) in chunk.enumerated() {
                    let base = 2 * p.grid.count * k
                    for i in 0..<p.grid.count {
                        v[base + 2 * i] = Float(f.real[i])
                        v[base + 2 * i + 1] = Float(f.imaginary[i])
                    }
                }
                let psi = try engine.floats(v, [chunk.count, p.grid.nx, p.grid.ny], .complexFloat32)
                try propagate(p, psi: psi, batch: chunk.count, planes: [], options: options, engine: engine) { j, data in
                    if j == p.sliceCount - 1 { out += fields(data, batch: chunk.count, grid: p.grid) }
                }
            }
        }
        return out
    }

    /// Slice by slice, V of a `Potential` as the GPU builds it, widened to float64 (for the tests of
    /// registration-performance.md section 3 and truth 1).
    package static func gpuPotentialSlices(_ potential: Potential) throws -> [[Double]] {
        let engine = try MetalEngine.require()
        engine.runLock.lock()
        defer { engine.runLock.unlock() }
        let beam = try ElectronBeam(energy: 100e3)
        let p = try prepare(potential, beam: beam, engine: engine)
        guard case .atoms(let e, let m) = p.source else { return [] }
        let nx = p.grid.nx, ny = p.grid.ny
        let step = engine.step("potential e\(e) m\(m) \(nx)x\(ny)") {
            let g = MPSGraph()
            let idx = g.placeholder(shape: [NSNumber(value: m)], dataType: .int32, name: "idx")
            let w = g.placeholder(shape: [NSNumber(value: m)], dataType: .float32, name: "w")
            let k = g.placeholder(shape: MetalGraphs.shape([e, nx, ny]), dataType: .float32, name: "kernels")
            let v = MetalGraphs.potential(g, idx: idx, w: w, kernels: k, elements: e, nx: nx, ny: ny)
            return MetalStep(executable: MetalExecutable(g, device: engine.graphDevice,
                                                         feeds: [(idx, [m], .int32), (w, [m], .float32),
                                                                 (k, [e, nx, ny], .float32)], targets: [v]),
                             placeholders: ["idx": idx, "w": w, "kernels": k])
        }
        return try (0..<p.sliceCount).map { j in
            let cb = MPSCommandBuffer(from: engine.queue)
            let out = step.executable.executable.encode(
                to: cb, inputs: step.inputs(["idx": p.indices[j], "w": p.weights[j], "kernels": p.kernels!]),
                results: nil, executionDescriptor: nil)
            cb.commit()
            cb.waitUntilCompleted()
            try engine.check(cb, "potential slice \(j)")
            var v = [Float](repeating: 0, count: p.grid.count)
            out[0].mpsndarray().readBytes(&v, strideBytes: nil)
            return v.map(Double.init)
        }
    }
}

extension MetalMultislice {
    /// IFFT(Ψ₀ · r_x ⊗ r_y) for a batch of b positions.
    static func probeStep(nx: Int, ny: Int, batch b: Int, engine: MetalEngine) -> MetalStep {
        return engine.step("probes \(nx)x\(ny) b\(b)") {
            let g = MPSGraph()
            let s = g.placeholder(shape: MetalGraphs.shape([nx, ny]), dataType: .complexFloat32, name: "spectrum")
            let rx = g.placeholder(shape: MetalGraphs.shape([b, nx, 1]), dataType: .complexFloat32, name: "rx")
            let ry = g.placeholder(shape: MetalGraphs.shape([b, 1, ny]), dataType: .complexFloat32, name: "ry")
            let shifted = g.multiplication(g.multiplication(s, rx, name: nil), ry, name: nil)
            let psi = g.fastFourierTransform(shifted, axes: [1, 2], descriptor: MetalGraphs.inverse, name: nil)
            return MetalStep(executable: MetalExecutable(g, device: engine.graphDevice, feeds: [
                (s, [nx, ny], .complexFloat32), (rx, [b, nx, 1], .complexFloat32), (ry, [b, 1, ny], .complexFloat32),
            ], targets: [psi]), placeholders: ["spectrum": s, "rx": rx, "ry": ry])
        }
    }

    /// The probes at the given ramps as the GPU builds them, widened to float64 (registration-performance.md 3).
    package static func gpuProbes(spectrum: ComplexField, ramps: [(x: ComplexField, y: ComplexField)],
                                  grid: Grid) throws -> [ComplexField] {
        let engine = try MetalEngine.require()
        engine.runLock.lock()
        defer { engine.runLock.unlock() }
        let b = ramps.count, nx = grid.nx, ny = grid.ny
        var rx = [Float](repeating: 0, count: 2 * b * nx), ry = [Float](repeating: 0, count: 2 * b * ny)
        for (k, r) in ramps.enumerated() {
            for i in 0..<nx { rx[2 * (k * nx + i)] = Float(r.x.real[i]); rx[2 * (k * nx + i) + 1] = Float(r.x.imaginary[i]) }
            for i in 0..<ny { ry[2 * (k * ny + i)] = Float(r.y.real[i]); ry[2 * (k * ny + i) + 1] = Float(r.y.imaginary[i]) }
        }
        let step = probeStep(nx: nx, ny: ny, batch: b, engine: engine)
        let (psi, _) = try engine.tensor([b, nx, ny], .complexFloat32)
        let cb = MPSCommandBuffer(from: engine.queue)
        _ = step.executable.executable.encode(
            to: cb, inputs: step.inputs(["spectrum": try engine.complex(spectrum),
                                         "rx": try engine.floats(rx, [b, nx, 1], .complexFloat32),
                                         "ry": try engine.floats(ry, [b, 1, ny], .complexFloat32)]),
            results: [psi], executionDescriptor: nil)
        cb.commit()
        cb.waitUntilCompleted()
        try engine.check(cb, "probes")
        return fields(psi, batch: b, grid: grid)
    }

    /// Scans (registration-performance.md 1.6): for each position the probe IFFT(Ψ₀ · r_x ⊗ r_y) built on the GPU,
    /// propagated, and its diffraction intensity |FFT(ψ)|² · scale read back widened to float64, `deliver(position, values)`.
    static func scanIntensities(spectrum: ComplexField, positions count: Int,
                                ramps: (Int) throws -> (x: ComplexField, y: ComplexField),
                                potential p: MetalPreparedPotential, intensityScale: Double,
                                        options: MetalOptions = .standard,
                                        deliver: (Int, [Double]) throws -> Void) throws {
        let engine = try MetalEngine.require()
        engine.runLock.lock()
        defer { engine.runLock.unlock() }
        let nx = p.grid.nx, ny = p.grid.ny, n = p.grid.count
        let size = batchSize(p.grid)
        let spectrumData = try engine.complex(spectrum)
        let scale = try engine.floats([Float(intensityScale)], [1])
        for start in stride(from: 0, to: count, by: size) {
            try autoreleasepool {
                let b = min(size, count - start)
                let probes = probeStep(nx: nx, ny: ny, batch: b, engine: engine)
                let intensities = engine.step("intensities \(nx)x\(ny) b\(b)") {
                    let g = MPSGraph()
                    let psi = g.placeholder(shape: MetalGraphs.shape([b, nx, ny]), dataType: .complexFloat32, name: "psi")
                    let sc = g.placeholder(shape: [1], dataType: .float32, name: "scale")
                    let f = g.fastFourierTransform(psi, axes: [1, 2], descriptor: MetalGraphs.forward, name: nil)
                    let re = g.realPartOfTensor(tensor: f, name: nil), im = g.imaginaryPartOfTensor(tensor: f, name: nil)
                    let i2 = g.addition(g.multiplication(re, re, name: nil), g.multiplication(im, im, name: nil), name: nil)
                    let out = g.multiplication(i2, sc, name: nil)
                    return MetalStep(executable: MetalExecutable(g, device: engine.graphDevice, feeds: [
                        (psi, [b, nx, ny], .complexFloat32), (sc, [1], .float32),
                    ], targets: [out]), placeholders: ["psi": psi, "scale": sc])
                }
                var rx = [Float](repeating: 0, count: 2 * b * nx), ry = [Float](repeating: 0, count: 2 * b * ny)
                for k in 0..<b {
                    let (fx, fy) = try ramps(start + k)
                    for i in 0..<nx {
                        rx[2 * (k * nx + i)] = Float(fx.real[i])
                        rx[2 * (k * nx + i) + 1] = Float(fx.imaginary[i])
                    }
                    for i in 0..<ny {
                        ry[2 * (k * ny + i)] = Float(fy.real[i])
                        ry[2 * (k * ny + i) + 1] = Float(fy.imaginary[i])
                    }
                }
                let rxData = try engine.floats(rx, [b, nx, 1], .complexFloat32)
                let ryData = try engine.floats(ry, [b, 1, ny], .complexFloat32)
                let (psi, _) = try engine.tensor([b, nx, ny], .complexFloat32)
                let cb = MPSCommandBuffer(from: engine.queue)
                _ = probes.executable.executable.encode(
                    to: cb, inputs: probes.inputs(["spectrum": spectrumData, "rx": rxData, "ry": ryData]), results: [psi],
                    executionDescriptor: nil)
                cb.commit()
                cb.waitUntilCompleted()
                try engine.check(cb, "probes")
                try propagate(p, psi: psi, batch: b, planes: [], options: options, engine: engine) { j, exit in
                    guard j == p.sliceCount - 1 else { return }
                    let (values, buffer) = try engine.tensor([b, nx, ny], .float32)
                    let cb2 = MPSCommandBuffer(from: engine.queue)
                    _ = intensities.executable.executable.encode(
                        to: cb2, inputs: intensities.inputs(["psi": exit, "scale": scale]), results: [values],
                        executionDescriptor: nil)
                    cb2.commit()
                    cb2.waitUntilCompleted()
                    try engine.check(cb2, "intensities")
                    let pointer = buffer.contents().bindMemory(to: Float.self, capacity: b * n)
                    for k in 0..<b {
                        try deliver(start + k, (0..<n).map { Double(pointer[k * n + $0]) })
                    }
                }
            }
        }
    }
}

extension Multislice {
    /// A scan on Metal (registration-performance.md 1.6), for the STEM driver: the probe Ψ₀ (`spectrum`, normalised,
    /// unshifted) shifted by `ramps(q)` and built on the GPU, propagated through `potential`, and the diffraction
    /// intensity |FFT(ψ)|² of each position delivered in position order, widened to float64.
    package static func metalScanIntensities(spectrum: ComplexField, positions: Int,
                                             ramps: (Int) throws -> (x: ComplexField, y: ComplexField),
                                             through potential: ProjectedPotential, beam: ElectronBeam,
                                             deliver: (Int, [Double]) throws -> Void) throws {
        let engine = try MetalEngine.require()
        let prepared = try MetalMultislice.prepare(potential, beam: beam, engine: engine)
        try MetalMultislice.scanIntensities(spectrum: spectrum, positions: positions, ramps: ramps, potential: prepared,
                                            intensityScale: 1, deliver: deliver)
    }

    /// The same through a `Potential` built on the GPU.
    package static func metalScanIntensities(spectrum: ComplexField, positions: Int,
                                             ramps: (Int) throws -> (x: ComplexField, y: ComplexField),
                                             through potential: Potential, beam: ElectronBeam,
                                             deliver: (Int, [Double]) throws -> Void) throws {
        let engine = try MetalEngine.require()
        let prepared = try MetalMultislice.prepare(potential, beam: beam, engine: engine)
        try MetalMultislice.scanIntensities(spectrum: spectrum, positions: positions, ramps: ramps, potential: prepared,
                                            intensityScale: 1, deliver: deliver)
    }
}
