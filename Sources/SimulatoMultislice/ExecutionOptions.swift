//
//  ExecutionOptions.swift
//  Role: the execution choices of both backends as parameters (registration-generality.md section 1): every default is
//        the value the kit was validated with, so `.standard` reproduces every registered number bit for bit; and the
//        registered agreement of `.metal` with `.cpu` as a public value (`MetalAccuracy.registered`).
//

/// How a remainder batch is shaped on Metal.
public enum MetalBatchShapes: Sendable, Hashable {
    /// Rounded up to a power of two (at most the batch size), the extra waves zero and not read: a run compiles a few
    /// shapes. The default.
    case powersOfTwo
    /// Exactly the remainder: one compiled shape per distinct batch size.
    case exact
}

/// How `Multislice.run`, `Multislice.runBatch` and `STEM.scan` execute. None of these choices changes a `.cpu` number;
/// on `.metal`, batch sizes and shapes change rounding only (within `MetalAccuracy.registered`), the renormalisation
/// changes the agreement (`docs/registration-performance.md`, amendment 1).
public struct ExecutionOptions: Sendable, Hashable {
    /// At most this many threads for a CPU batch or scan (default 8); each thread runs whole positions, bit for bit as a
    /// sequential run.
    public var maximumCPUThreads: Int
    /// The bytes of waves one Metal batch may hold (default 512 MiB; a wave costs nx · ny · 24 bytes).
    public var metalBatchBytes: Int
    /// The bytes of exit waves a CPU scan holds at once (default 256 MiB; an exit wave costs nx · ny · 32 bytes).
    public var scanChunkBytes: Int
    /// Rescale each Metal step to the device-measured Parseval energies (default on).
    public var metalRenormalization: Bool
    /// The compiled Metal executables kept, least recently used dropped first (default 12; each holds device memory).
    public var metalExecutableCacheLimit: Int
    /// How a remainder batch is shaped on Metal (default `.powersOfTwo`).
    public var metalBatchShapes: MetalBatchShapes

    public init(maximumCPUThreads: Int = 8, metalBatchBytes: Int = 512 << 20, scanChunkBytes: Int = 256 << 20,
                metalRenormalization: Bool = true, metalExecutableCacheLimit: Int = 12,
                metalBatchShapes: MetalBatchShapes = .powersOfTwo) {
        self.maximumCPUThreads = max(1, maximumCPUThreads)
        self.metalBatchBytes = max(1, metalBatchBytes)
        self.scanChunkBytes = max(1, scanChunkBytes)
        self.metalRenormalization = metalRenormalization
        self.metalExecutableCacheLimit = max(1, metalExecutableCacheLimit)
        self.metalBatchShapes = metalBatchShapes
    }

    /// The defaults the kit was validated with.
    public static let standard = ExecutionOptions()
}

/// The agreement of `.metal` (float32) with `.cpu` (float64) the kit registers and tests, relative to the array's
/// largest magnitude (`docs/registration-slice-1.md` section 4, amendment 2; `docs/registration-slice-2.md` section 4;
/// `docs/registration-performance.md`). Signals far below the array's maximum, under about 1e−5 of the incident
/// intensity, sit at float32's floor and are not covered (performance amendment 3, item 6).
public struct MetalAccuracy: Sendable, Hashable {
    /// Exit waves, ×max |ψ|.
    public var wave: Double
    /// Intensities |ψ|², ×max.
    public var intensity: Double
    /// Diffraction patterns, ×max.
    public var pattern: Double
    /// Detector values over a scan, ×the image's max.
    public var detector: Double
    /// Σ|ψ|², relative, per slice.
    public var totalIntensityPerSlice: Double
    /// The bounds hold up to this many slices.
    public var sliceLimit: Int
    /// And on grids up to this many points per side.
    public var gridLimit: Int

    public init(wave: Double = 1e-4, intensity: Double = 2e-4, pattern: Double = 2e-4, detector: Double = 2e-4,
                totalIntensityPerSlice: Double = 1.5e-6, sliceLimit: Int = 100, gridLimit: Int = 1024) {
        self.wave = wave
        self.intensity = intensity
        self.pattern = pattern
        self.detector = detector
        self.totalIntensityPerSlice = totalIntensityPerSlice
        self.sliceLimit = sliceLimit
        self.gridLimit = gridLimit
    }

    /// The registered bounds.
    public static let registered = MetalAccuracy()
}
