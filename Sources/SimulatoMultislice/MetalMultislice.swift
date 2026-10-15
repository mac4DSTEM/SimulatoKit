//
//  MetalMultislice.swift
//  Role: the backends' entry points (registration section 4, second table; registration-performance.md): `.metal`
//        runs on the GPU through `MetalEngine.swift` (float32, the potential and transmissions on the device); `.cpu`
//        batches run positions on at most 8 threads, bit for bit equal to sequential runs.
//

import Foundation
import Metal
import Synchronization
import SimulatoCore

enum MetalMultislice {
    /// Whether a Metal device exists (the tests skip the Metal path with a named reason when it does not).
    static var isAvailable: Bool { MetalEngine.shared != nil }
}

extension Multislice {
    /// Whether `.metal` can run here.
    public static var metalIsAvailable: Bool { MetalMultislice.isAvailable }
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

    /// The exit waves (after the last slice) of a batch of incident waves through one potential: on `.cpu` on at most
    /// 8 threads from the shared transmission functions, on `.metal` in batches whose waves fit in 512 MB
    /// (registration-slice-2.md 3.6), the transmissions built on the GPU.
    public static func runBatch(_ incidents: [Waves], through potential: ProjectedPotential,
                                backend: Backend = .cpu) throws -> [Waves] {
        guard let first = incidents.first else { return [] }
        for w in incidents where !w.grid.matches(potential.grid) { throw SimulatoError.gridMismatch }
        switch backend {
        case .cpu:
            let transmissions = try transmissionFunctions(potential, beam: first.beam)
            return try runCPUParallel(incidents, potential: potential, transmissions: transmissions,
                                      last: potential.slices.count - 1)
        case .metal:
            let engine = try MetalEngine.require()
            let prepared = try MetalMultislice.prepare(potential, beam: first.beam, engine: engine)
            return try metalBatch(incidents, prepared)
        }
    }

    /// The same for a `Potential`: on `.cpu` built first (`build()`), on `.metal` built slice by slice on the GPU.
    public static func runBatch(_ incidents: [Waves], through potential: Potential,
                                backend: Backend = .cpu) throws -> [Waves] {
        guard let first = incidents.first else { return [] }
        for w in incidents where !w.grid.matches(potential.grid) { throw SimulatoError.gridMismatch }
        switch backend {
        case .cpu:
            return try runBatch(incidents, through: potential.build(), backend: .cpu)
        case .metal:
            let engine = try MetalEngine.require()
            let prepared = try MetalMultislice.prepare(potential, beam: first.beam, engine: engine)
            return try metalBatch(incidents, prepared)
        }
    }

    static func metalBatch(_ incidents: [Waves], _ prepared: MetalPreparedPotential) throws -> [Waves] {
        let fields = try MetalMultislice.runBatch(incidents.map(\.field), prepared: prepared)
        return zip(incidents, fields).map {
            Waves(field: $0.1, grid: prepared.grid, beam: $0.0.beam, normalization: $0.0.normalization)
        }
    }
}

extension Multislice {
    /// For the tests of registration-performance.md: V of every slice of `potential` as the GPU builds it.
    package static func gpuPotentialSlices(_ potential: Potential) throws -> [[Double]] {
        try MetalMultislice.gpuPotentialSlices(potential)
    }

    /// For the tests: the probes Ψ₀ shifted by the ramps, built on the GPU.
    package static func gpuProbes(spectrum: ComplexField, ramps: [(x: ComplexField, y: ComplexField)],
                                  grid: Grid) throws -> [ComplexField] {
        try MetalMultislice.gpuProbes(spectrum: spectrum, ramps: ramps, grid: grid)
    }

    /// For the tests: `.metal` with the per-step renormalisation on or off (registration-performance.md, amendment 1).
    package static func runMetal(_ incident: Waves, through potential: Potential, exitPlanes: ExitPlanes = .last,
                                 renormalize: Bool) throws -> [Waves] {
        let engine = try MetalEngine.require()
        let prepared = try MetalMultislice.prepare(potential, beam: incident.beam, engine: engine)
        return try MetalMultislice.run(incident, prepared: prepared,
                                       planes: exitPlanes.planes(sliceCount: potential.sliceCount),
                                       options: MetalOptions(renormalize: renormalize))
    }

    package static func runMetal(_ incident: Waves, through potential: ProjectedPotential, exitPlanes: ExitPlanes = .last,
                                 renormalize: Bool) throws -> [Waves] {
        let engine = try MetalEngine.require()
        let prepared = try MetalMultislice.prepare(potential, beam: incident.beam, engine: engine)
        return try MetalMultislice.run(incident, prepared: prepared,
                                       planes: exitPlanes.planes(sliceCount: potential.slices.count),
                                       options: MetalOptions(renormalize: renormalize))
    }
}

extension Multislice {
    /// The Metal device's largest allocated size [bytes] seen at a commit since the last reset (for the memory report of
    /// registration-performance.md 1.8); `reset` starts a new measurement at the current size.
    package static func metalPeakAllocatedBytes(reset: Bool = false) -> Int {
        guard let engine = MetalEngine.shared else { return 0 }
        if reset { engine.peak.withLock { $0 = engine.device.currentAllocatedSize } }
        return engine.peak.withLock { $0 }
    }
}
