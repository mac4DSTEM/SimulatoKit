//
//  MetalMultislice.swift
//  Role: the accelerated path (registration section 4, second table): the same multislice loop in float32 on the GPU,
//        one MPSGraph per grid: out = IFFT(FFT(ψ · t) · P) with the inverse scaled by 1/N (NumPy's convention). The
//        bandlimited transmission functions and the propagators are computed by the CPU path in float64 and cast to
//        complex float32 once; the exit waves are read back and widened to float64.
//

import Foundation
import Metal
import Synchronization
import MetalPerformanceShadersGraph
import SimulatoCore

enum MetalMultislice {
    /// Whether a Metal device exists (the tests skip the Metal path with a named reason when it does not).
    static var isAvailable: Bool { MTLCreateSystemDefaultDevice() != nil }

    static func run(_ incident: Waves, potential: ProjectedPotential, transmissions: [ComplexField],
                    planes: [Int]) throws -> [Waves] {
        guard let device = MTLCreateSystemDefaultDevice() else {
            throw SimulatoError.metalUnavailable("no Metal device (MTLCreateSystemDefaultDevice() returned nil)")
        }
        guard let queue = device.makeCommandQueue() else {
            throw SimulatoError.metalUnavailable("the device made no command queue")
        }
        let grid = potential.grid
        let shape: [NSNumber] = [NSNumber(value: grid.nx), NSNumber(value: grid.ny)]
        let graph = MPSGraph()
        let psiIn = graph.placeholder(shape: shape, dataType: .complexFloat32, name: "psi")
        let tIn = graph.placeholder(shape: shape, dataType: .complexFloat32, name: "t")
        let pIn = graph.placeholder(shape: shape, dataType: .complexFloat32, name: "p")
        let forward = MPSGraphFFTDescriptor()
        forward.inverse = false
        forward.scalingMode = .none
        let inverse = MPSGraphFFTDescriptor()
        inverse.inverse = true
        inverse.scalingMode = .size
        let transmitted = graph.multiplication(psiIn, tIn, name: nil)
        let spectrum = graph.fastFourierTransform(transmitted, axes: [0, 1], descriptor: forward, name: nil)
        let propagated = graph.multiplication(spectrum, pIn, name: nil)
        let out = graph.fastFourierTransform(propagated, axes: [0, 1], descriptor: inverse, name: nil)

        let graphDevice = MPSGraphDevice(mtlDevice: device)
        func tensorData(_ f: ComplexField) -> MPSGraphTensorData {
            var interleaved = [Float](repeating: 0, count: 2 * f.real.count)
            for i in f.real.indices {
                interleaved[2 * i] = Float(f.real[i])
                interleaved[2 * i + 1] = Float(f.imaginary[i])
            }
            return MPSGraphTensorData(device: graphDevice, data: Data(bytes: interleaved, count: interleaved.count * 4),
                                      shape: shape, dataType: .complexFloat32)
        }
        func field(_ data: MPSGraphTensorData) -> ComplexField {
            var interleaved = [Float](repeating: 0, count: 2 * grid.count)
            data.mpsndarray().readBytes(&interleaved, strideBytes: nil)
            var f = ComplexField(nx: grid.nx, ny: grid.ny)
            for i in 0..<grid.count {
                f.real[i] = Double(interleaved[2 * i])
                f.imaginary[i] = Double(interleaved[2 * i + 1])
            }
            return f
        }

        let tData = transmissions.map(tensorData)
        var pData: [Double: MPSGraphTensorData] = [:]
        var psi = tensorData(incident.field)
        var result: [Waves] = []
        if planes.first == -1 { result.append(incident) }
        let wanted = Set(planes)
        for (j, slice) in potential.slices.enumerated() {
            if pData[slice.thickness] == nil {
                pData[slice.thickness] = tensorData(Multislice.propagator(grid: grid, beam: incident.beam,
                                                                          thickness: slice.thickness))
            }
            let feeds: [MPSGraphTensor: MPSGraphTensorData] = [psiIn: psi, tIn: tData[j], pIn: pData[slice.thickness]!]
            guard let next = graph.run(with: queue, feeds: feeds, targetTensors: [out], targetOperations: nil)[out] else {
                throw SimulatoError.metalFailure("the graph returned no output at slice \(j)")
            }
            psi = next
            if wanted.contains(j) {
                result.append(Waves(field: field(psi), grid: grid, beam: incident.beam,
                                    normalization: incident.normalization))
            }
        }
        return result
    }
}

extension Multislice {
    /// Whether `.metal` can run here.
    public static var metalIsAvailable: Bool { MetalMultislice.isAvailable }
}

extension MetalMultislice {
    /// Slice 2 (registration-slice-2.md 3.6): a batch of incident waves through one potential, ψ of shape
    /// [B, nx, ny] complex64, t and P broadcast, the FFT over axes 1 and 2; the exit waves after the last slice, in input
    /// order, widened to float64.
    package static func runBatch(_ incidents: [ComplexField], potential: ProjectedPotential, beam: ElectronBeam,
                                 transmissions: [ComplexField]) throws -> [ComplexField] {
        guard !incidents.isEmpty else { return [] }
        guard let device = MTLCreateSystemDefaultDevice() else {
            throw SimulatoError.metalUnavailable("no Metal device (MTLCreateSystemDefaultDevice() returned nil)")
        }
        guard let queue = device.makeCommandQueue() else {
            throw SimulatoError.metalUnavailable("the device made no command queue")
        }
        let grid = potential.grid
        let b = incidents.count
        let waveShape: [NSNumber] = [NSNumber(value: b), NSNumber(value: grid.nx), NSNumber(value: grid.ny)]
        let planeShape: [NSNumber] = [NSNumber(value: grid.nx), NSNumber(value: grid.ny)]
        let graph = MPSGraph()
        let psiIn = graph.placeholder(shape: waveShape, dataType: .complexFloat32, name: "psi")
        let tIn = graph.placeholder(shape: planeShape, dataType: .complexFloat32, name: "t")
        let pIn = graph.placeholder(shape: planeShape, dataType: .complexFloat32, name: "p")
        let forward = MPSGraphFFTDescriptor()
        forward.inverse = false
        forward.scalingMode = .none
        let inverse = MPSGraphFFTDescriptor()
        inverse.inverse = true
        inverse.scalingMode = .size
        let transmitted = graph.multiplication(psiIn, tIn, name: nil)
        let spectrum = graph.fastFourierTransform(transmitted, axes: [1, 2], descriptor: forward, name: nil)
        let propagated = graph.multiplication(spectrum, pIn, name: nil)
        let out = graph.fastFourierTransform(propagated, axes: [1, 2], descriptor: inverse, name: nil)
        let graphDevice = MPSGraphDevice(mtlDevice: device)
        func data(_ fields: [ComplexField], shape: [NSNumber]) -> MPSGraphTensorData {
            var interleaved = [Float](repeating: 0, count: 2 * grid.count * fields.count)
            for (k, f) in fields.enumerated() {
                let base = 2 * grid.count * k
                for i in 0..<grid.count {
                    interleaved[base + 2 * i] = Float(f.real[i])
                    interleaved[base + 2 * i + 1] = Float(f.imaginary[i])
                }
            }
            return MPSGraphTensorData(device: graphDevice, data: Data(bytes: interleaved, count: interleaved.count * 4),
                                      shape: shape, dataType: .complexFloat32)
        }
        var psi = data(incidents, shape: waveShape)
        var pData: [Double: MPSGraphTensorData] = [:]
        for (j, slice) in potential.slices.enumerated() {
            if pData[slice.thickness] == nil {
                pData[slice.thickness] = data([Multislice.propagator(grid: grid, beam: beam, thickness: slice.thickness)],
                                              shape: planeShape)
            }
            let feeds: [MPSGraphTensor: MPSGraphTensorData] = [psiIn: psi, tIn: data([transmissions[j]], shape: planeShape),
                                                               pIn: pData[slice.thickness]!]
            guard let next = graph.run(with: queue, feeds: feeds, targetTensors: [out], targetOperations: nil)[out] else {
                throw SimulatoError.metalFailure("the batched graph returned no output at slice \(j)")
            }
            psi = next
        }
        var interleaved = [Float](repeating: 0, count: 2 * grid.count * b)
        psi.mpsndarray().readBytes(&interleaved, strideBytes: nil)
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
}

extension Multislice {
    /// The CPU batch on at most `maximumCPUThreads` threads (the suite's 8-job cap): each thread runs whole positions
    /// with its own FFT plan; every exit wave is computed exactly as one sequential run computes it (bit for bit).
    public static let maximumCPUThreads = 8

    static func runCPUParallel(_ incidents: [Waves], potential: ProjectedPotential, transmissions: [ComplexField],
                               last: Int) throws -> [Waves] {
        let count = incidents.count
        let workers = max(1, min(maximumCPUThreads, ProcessInfo.processInfo.activeProcessorCount, count))
        let results = Mutex<[Int: Waves]>([:])
        let failure = Mutex<(any Error)?>(nil)
        DispatchQueue.concurrentPerform(iterations: workers) { worker in
            for q in stride(from: worker, to: count, by: workers) {
                do {
                    let w = try runCPU(incidents[q], potential: potential, transmissions: transmissions, planes: [last])[0]
                    results.withLock { $0[q] = w }
                } catch {
                    failure.withLock { if $0 == nil { $0 = error } }
                    return
                }
            }
        }
        if let error = failure.withLock({ $0 }) { throw error }
        return results.withLock { r in (0..<count).map { r[$0]! } }
    }

    /// The exit waves (after the last slice) of a batch of incident waves through one potential: on `.cpu` one at a time
    /// from the shared transmission functions, on `.metal` in batches whose waves fit in 512 MB (registration-slice-2.md
    /// 3.6). The transmission functions are computed once.
    public static func runBatch(_ incidents: [Waves], through potential: ProjectedPotential,
                                backend: Backend = .cpu) throws -> [Waves] {
        guard let first = incidents.first else { return [] }
        for w in incidents where !w.grid.matches(potential.grid) { throw SimulatoError.gridMismatch }
        let transmissions = try transmissionFunctions(potential, beam: first.beam)
        let last = potential.slices.count - 1
        switch backend {
        case .cpu:
            return try runCPUParallel(incidents, potential: potential, transmissions: transmissions, last: last)
        case .metal:
            let perWave = potential.grid.count * 8 * 3  // the wave, its spectrum and the product, complex64
            let batch = max(1, (512 << 20) / perWave)
            var out: [Waves] = []
            for start in stride(from: 0, to: incidents.count, by: batch) {
                let chunk = Array(incidents[start..<min(start + batch, incidents.count)])
                let fields = try MetalMultislice.runBatch(chunk.map(\.field), potential: potential, beam: first.beam,
                                                          transmissions: transmissions)
                out += zip(chunk, fields).map {
                    Waves(field: $0.1, grid: potential.grid, beam: $0.0.beam, normalization: $0.0.normalization)
                }
            }
            return out
        }
    }
}
