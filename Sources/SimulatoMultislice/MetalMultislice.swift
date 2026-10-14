//
//  MetalMultislice.swift
//  Role: the accelerated path (registration section 4, second table): the same multislice loop in float32 on the GPU,
//        one MPSGraph per grid: out = IFFT(FFT(ψ · t) · P) with the inverse scaled by 1/N (NumPy's convention). The
//        bandlimited transmission functions and the propagators are computed by the CPU path in float64 and cast to
//        complex float32 once; the exit waves are read back and widened to float64.
//

import Foundation
import Metal
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
