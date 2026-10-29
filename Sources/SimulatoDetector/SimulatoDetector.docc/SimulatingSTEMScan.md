# Simulating a STEM scan

Scan a probe over a specimen, record several detectors at once, average over frozen-phonon configurations, and write a detector of your own.

## Overview

`STEM.scan` is the one entry point. It takes a `Probe`, the potentials to scan through, a scan (`GridScan` or `LineScan`) and a list of detectors, and returns one ``ScanMeasurement`` per detector, in the order given. A scan computes every exit wave once and feeds all detectors from it.

### A HAADF image

```swift
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe
import SimulatoDetector

let beam = try ElectronBeam(energy: 300e3)
let potential = try Potential(structure: structure, nx: 256, ny: 256, sliceThickness: 2.0)
let probe = try Probe(beam: beam, grid: potential.grid,
                      aperture: Aperture(semiangleCutoff: 20))
let scan = try GridScan(start: SIMD2(0, 0),
                        end: SIMD2(structure.box.x, structure.box.y), nx: 32, ny: 32)

let backend: Backend = Multislice.metalIsAvailable ? .metal : .cpu
let measurements = try STEM.scan(probe: probe, potentials: [potential], scan: scan,
                                 detectors: [AnnularDetector(inner: 70, outer: 200)],
                                 backend: backend)
let haadf = measurements[0].image(samplingX: structure.box.x / 32,
                                  samplingY: structure.box.y / 32)!
```

The probe and every potential must be on the same grid (`Grid.matches(_:)`); otherwise `STEM.scan` throws `SimulatoError.gridMismatch`. The potentials can be built (`ProjectedPotential`) or unbuilt (`Potential`). With unbuilt potentials `.cpu` builds each one first and `.metal` builds each slice on the GPU.

### Detectors

All angles are in mrad.

- ``AnnularDetector``: sums the diffraction intensity where inner² ≤ α² < outer². One number per position. `outer` defaults to the largest whole angle inside the antialias cutoff, and an outer angle beyond the simulated pattern throws. A bright-field detector is `AnnularDetector(inner: 0, outer: semiangle)`, and ADF and HAADF differ only in the limits.
- ``PixelatedDetector``: the diffraction pattern, cropped and shifted so that zero angle is at the centre. `maxAngle` is `.cutoff`, `.valid` (the default, the region free of aliasing, which is the cutoff divided by √2) or `.full`. The value shape is `[nx, ny]` of the cropped pattern.
- ``SegmentedDetector``: `radialBins` x `azimuthalBins` bins between `inner` and `outer`, the azimuth measured from `rotation` in rad. The value shape is `[radialBins, azimuthalBins]`, with the azimuthal index fastest.
- ``CenterOfMassDetector``: the centre of mass of the pattern over a pixelated detector's region, in 1/Å (default) or mrad. Two values per position, x then y.
- ``WavesDetector``: the complex exit wave itself, optionally resampled in Fourier space to `rows` x `columns`. It has its own `STEM.scan` overload and returns ``ScanWaves`` (see below).

### ScanMeasurement layout

``ScanMeasurement`` stores `values` row-major over (scan x, scan y, the detector's value shape), with the detector's own values contiguous at each position:

```
values[(ix * scanShape[1] + iy) * valueCount + k]
```

`scanShape` is `[scan.nx, scan.ny]` (or `[n]` for a `LineScan`), `valueShape` is the detector's shape (empty for a single number), and `subscript(ix, iy)` returns the values at one position. `image(samplingX:samplingY:)` wraps a one-number detector's values as a `RealField` in the scan's x-first layout and returns `nil` for a detector with more values.

```swift
let pixelated = try STEM.scan(probe: probe, potentials: [potential], scan: scan,
                              detectors: [PixelatedDetector(maxAngle: .valid)])[0]
let pattern = pixelated[3, 5]                    // one diffraction pattern
let (px, py) = (pixelated.valueShape[0], pixelated.valueShape[1])
```

### Line scans and profiles

A `LineScan` places `gpts` positions (or one every `sampling` Å) on a straight line. The measurement has one scan axis (`scanShape == [n]`), and `lineProfiles(for:)` turns it into ``RealSpaceLineProfiles``, with `width(height:)` (the full width at a fraction of the maximum) and `interpolate(sampling:gpts:order:endpoint:)`.

```swift
let line = try LineScan(start: SIMD2(0, 0), end: SIMD2(structure.box.x, 0), gpts: 64)
let adf = try STEM.scan(probe: probe, potentials: [potential], scan: line,
                        detectors: [AnnularDetector(inner: 50, outer: 200)])[0]
let profile = try adf.lineProfiles(for: line)
let fwhm = try profile.width()                    // one width per profile, in Å
```

A profile of an image you already have comes from `RealField.interpolateLine(start:end:sampling:gpts:...)`, which samples along the line with spline orders 0 to 5, as abTEM's `interpolate_line`.

### Integrated centre of mass

A ``CenterOfMassDetector`` scan gives the gradient field; `integratedCenterOfMass(samplingX:samplingY:)` integrates it in Fourier space (abTEM's `integrated_center_of_mass`). `gradientImages(samplingX:samplingY:)` returns the two components as images.

```swift
let com = try STEM.scan(probe: probe, potentials: [potential], scan: scan,
                        detectors: [CenterOfMassDetector()])[0]
let icom = try com.integratedCenterOfMass(samplingX: structure.box.x / 32, samplingY: structure.box.y / 32)
```

### Exit waves at every position

`STEM.scan(probe:potentials:scan:detector:backend:options:memoryLimit:)` with a ``WavesDetector`` returns ``ScanWaves``: the complex exit wave at every position of a `GridScan`, `scanWaves[ix, iy]` as a `Waves`. It takes one potential (frozen phonons are refused with `WavesDetectorError.severalPotentials`; scan each configuration on its own) and refuses a result larger than half of the physical memory or the `memoryLimit` given (`WavesDetectorError.tooLarge`).

### Frozen-phonon averaging

Thermal displacement is simulated by averaging intensities over configurations of the structure. `FrozenPhonons` draws a displacement for every atom from a Gaussian of standard deviation σ per element (in Å), with a seeded generator, so the same seed gives the same ensemble on every platform. Pass the resulting potentials to `STEM.scan` and each detector value is averaged over them.

```swift
import SimulatoMultislice

let phonons = try FrozenPhonons(structure: structure,
                                sigmas: [13: 0.08], configurations: 8, seed: 1)
let potentials = try phonons.lazyPotentials(nx: 256, ny: 256, sliceThickness: 2.0)
let averaged = try STEM.scan(probe: probe, potentials: potentials, scan: scan,
                             detectors: [AnnularDetector(inner: 70, outer: 200)],
                             backend: .metal)[0]
```

(`sigmas` maps the atomic number to σ; every element of the structure needs an entry. Use `potentials(nx:ny:sliceThickness:)` instead for built potentials.) The average is of intensities, taken after the detector: each position's values are summed over configurations in potential order and divided by the number of potentials.

### Writing your own detector

A custom detector conforms to ``Detector``: `shape(for:beam:)` gives the shape of one position's values (`[]` for a single number), and `detect(_:)` reduces a ``DiffractionIntensity`` to that many values, row-major. `STEM.scan` checks that the returned count equals the shape's product.

``DiffractionIntensity`` gives the full unshifted pattern |FFT(ψ)|² in FFT order, `values[ix * ny + iy]`, together with the grid and beam. Its `angularSampling` (mrad per pixel) and `Grid.fftfreq(_:_:)` convert indices to angles, as in the example.

```swift
/// Sum over a disc offset from the optic axis (for example one DPC-style segment).
struct OffsetDisc: Detector {
    let centre: SIMD2<Double>      // mrad
    let radius: Double             // mrad

    func shape(for grid: Grid, beam: ElectronBeam) throws -> [Int] { [] }

    func detect(_ p: DiffractionIntensity) throws -> [Double] {
        let a = p.angularSampling
        let kx = Grid.fftfreq(p.grid.nx, 1 / a.x / Double(p.grid.nx))
        let ky = Grid.fftfreq(p.grid.ny, 1 / a.y / Double(p.grid.ny))
        var sum = 0.0
        for ix in 0..<p.grid.nx {
            for iy in 0..<p.grid.ny {
                let dx = kx[ix] - centre.x, dy = ky[iy] - centre.y
                if dx * dx + dy * dy < radius * radius { sum += p.values[ix * p.grid.ny + iy] }
            }
        }
        return [sum]
    }
}
```

Test a detector without running a simulation by constructing the intensity directly:

```swift
let grid = try Grid(nx: 64, ny: 64, extentX: 10, extentY: 10)
let beam = try ElectronBeam(energy: 300e3)
let flat = DiffractionIntensity(grid: grid, beam: beam,
                                values: [Double](repeating: 1, count: grid.count))
let value = try OffsetDisc(centre: .zero, radius: 10).detect(flat)
```

### Backends and memory

`STEM.scan` also takes an `options:` parameter of type `ExecutionOptions` (SimulatoMultislice). It sets the CPU thread cap, the Metal batch memory and the size of the chunk of positions held at once (`scanChunkBytes`); `ExecutionOptions.standard` is the default.

On `.cpu` the positions run in chunks and in parallel on at most 8 threads. On `.metal` the probes, the transmissions and the diffraction intensities are computed on the GPU in batches (see `Metal on Apple GPUs: pitfalls and design` in SimulatoMultislice). Results from the two backends agree within the registered bounds; use `.cpu` as the reference in tests.
