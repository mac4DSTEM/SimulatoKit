# SimulatoKit

A Swift package for simulating electron scattering on Apple platforms: projected potentials, multislice propagation,
frozen phonons, probes with aberrations, STEM scans with detectors, and the contrast transfer function for HRTEM images.
It is a port of [abTEM](https://github.com/abtem/abtem) 1.0.10, checked against abTEM in float64 (Python oracle, pinned
inputs, `docs/registration-slice-*.md`). The algorithms and the Lobato–Van Dyck (2014) parameters are abTEM's; credit
for the method belongs to the abTEM authors. Licence: GPL-3.0-or-later (`LICENSE`, `NOTICE`); software built on
SimulatoKit must be GPL-compatible. It is a library of numbers, not an application, and has no UI.

## Modules

| Product | What it does |
|---|---|
| `SimulatoCore` | the electron beam (wavelength, interaction parameter), grids, complex and real fields, waves (intensity, diffraction pattern), the antialias aperture, a minimal orthogonal-cell `AtomicStructure`, a float64 FFT |
| `SimulatoMultislice` | the Lobato parameterisation, slicing, the projected potential, transmission and Fresnel propagation, `Multislice.run` / `runBatch` on `.cpu` and `.metal`, frozen phonons |
| `SimulatoProbe` | the soft aperture, polar aberrations C10 to C56, the probe at a position, `GridScan`, `ContrastTransferFunction` (envelopes, Scherzer defocus, point resolution) and `Waves.applying` for HRTEM |
| `SimulatoDetector` | annular (BF, ADF, HAADF), pixelated (4D-STEM) and segmented (DPC) detectors, and `STEM.scan`, averaged over frozen-phonon configurations |

## Installation

```swift
.package(url: "https://github.com/mac4DSTEM/SimulatoKit", from: "0.1.0")
// target dependencies, as needed:
.product(name: "SimulatoMultislice", package: "SimulatoKit")  // also SimulatoCore, SimulatoProbe, SimulatoDetector
```

Swift 6.2 tools, Swift 6 language mode. Platforms: macOS 15 and iOS/iPadOS 18 (from 0.2.0). Every module, the Metal
path included, builds for both (MPSGraph's FFT needs macOS 14 / iOS 17, `Synchronization.Mutex` macOS 15 / iOS 18); the
tests and the consumer checks run on macOS, where the numbers are validated. Developed and tested on macOS 27.
Build with `-c release`: the FFT and potential loops are many times slower unoptimised. Grid lengths are f · 2^m with
f ∈ {1, 3, 5, 15} and m ≥ 3 (the lengths Accelerate's DFT accepts). Cells are orthogonal.

## Quick start

```swift
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe
import SimulatoDetector

let a = 4.05, h = a / 2   // fcc Al, Å
let cell = try AtomicStructure(box: SIMD3(a, a, a), atoms: [
    .init(atomicNumber: 13, position: SIMD3(0, 0, 0)), .init(atomicNumber: 13, position: SIMD3(0, h, h)),
    .init(atomicNumber: 13, position: SIMD3(h, 0, h)), .init(atomicNumber: 13, position: SIMD3(h, h, 0))])
let beam = try ElectronBeam(energy: 300e3)   // eV
let backend: Backend = Multislice.metalIsAvailable ? .metal : .cpu
// A Potential is the sliced structure: built on the GPU slice by slice on .metal, or on the CPU (`build()`) on .cpu.
let potential = try Potential(structure: cell.repeated(3, 3, 10), nx: 256, ny: 256, sliceThickness: 2)

// HRTEM: plane wave, multislice, CTF at Scherzer defocus (Cs = 1 mm; lengths in Å)
let exit = try Multislice.run(.planeWave(beam: beam, grid: potential.grid), through: potential, backend: backend)[0]
var ab = Aberrations()
ab.Cs = 1e7
ab.defocus = try ContrastTransferFunction.scherzerDefocus(Cs: 1e7, energy: 300e3)
let hrtem = try exit.applying(ContrastTransferFunction(aberrations: ab)).intensity()   // RealField

// STEM HAADF: 20 mrad probe over one unit cell, annular detector 70 to 200 mrad
let probe = try Probe(beam: beam, grid: potential.grid, aperture: Aperture(semiangleCutoff: 20))
let scan = try GridScan(start: SIMD2(0, 0), end: SIMD2(a, a), nx: 16, ny: 16)
let haadf = try STEM.scan(probe: probe, potentials: [potential], scan: scan,
                          detectors: [AnnularDetector(inner: 70, outer: 200)],
                          backend: backend)[0]   // ScanMeasurement: values over the scan
```

## Backends

`Backend.cpu` is the reference: float64 on Accelerate (vDSP); scans and batches use up to 8 threads. `Backend.metal`
runs in float32 on the GPU (MPSGraph) behind the same API: for a `Potential` the projected potential is built on the GPU
slice by slice, the transmission functions always are, the wave stays on the device, and scans build their probes and
diffraction intensities there too; one compiled step is encoded per slice into one command buffer. Without a Metal
device it throws `SimulatoError.metalUnavailable` (`Multislice.metalIsAvailable` tells beforehand). The registered
tolerances of Metal against the CPU: exit wave 1e−4, intensity, diffraction pattern and detector images 2e−4 relative to
the maximum, total intensity 1.5e−6 per slice, up to 100 slices (`docs/registration-performance.md`); measured, the wave
is within about 5e−5 and the total intensity within 5e−7 per slice (1.4e−7 on light elements) at up to 1024² with 100
slices. Very weak high-angle signals (below about 1e−5 of the probe's intensity) sit at float32's noise floor.

Measured speed (Apple M5 Pro, macOS 27; fcc Al at 300 keV, 100 slices, about 0.047 Å per pixel; from the structure to the
numbers, the potential included; `docs/registration-performance.md`):

| Grid | Plane wave exit: Metal / CPU | 64 × 64 HAADF scan: Metal / CPU | Metal peak memory |
|---|---|---|---|
| 256² | 0.09 s / 0.59 s | 17 s / 516 s | 1.0 GB |
| 512² | 0.13 s / 4.4 s | 90 s / 2521 s | 1.3 GB |
| 1024² | 0.28 s / 54 s | 408 s / 15157 s | 1.7 GB |

The CPU scans run 8 threads and were timed on a subset of positions and scaled (linear in positions). These first
numbers were taken on a heavily loaded machine (other builds running), which slows the CPU side most: read the ratios
as upper bounds until a quiet re-measure.

## Choices and their defaults

Every choice the kit makes has a parameter, and the default is the value the numbers were validated with:
`ExecutionOptions` (`options:` on `Multislice.run`, `Multislice.runBatch` and `STEM.scan`: CPU threads, Metal batch memory
and shapes, the CPU scan's chunk, the Metal renormalisation and executable cache), `MetalAccuracy.registered` (the
Metal-against-CPU bounds the tests hold), `Grid.matches(_:relativeTolerance:)`, the slice ceiling (`maximumSliceCount:`)
and the constants of `scherzerDefocus`, `pointResolution` and `nyquistSampling`. abTEM's own conventions (the antialias
aperture, the slicing rules, the π choice in `LobatoParameterization`) stay as abTEM defines them.

```swift
// continuing the quick start: its 256 probe positions as a batch, with smaller Metal batches and fewer CPU threads
var options = ExecutionOptions.standard
options.maximumCPUThreads = 4
options.metalBatchBytes = 128 << 20
let probes = try scan.positions.map { try probe.build(at: $0) }
let exits = try Multislice.runBatch(probes, through: potential, backend: backend, options: options)
```

## How the numbers are checked

- **Pins against abTEM 1.0.10** (`164e644`) in float64: wavelength, potentials, probes, exit waves, patterns, detectors,
  CTF kernels and HRTEM images agree within 1e−10 of the array's maximum or tighter (measured 1e−13 to 1e−16, with
  abTEM's float32 π reproduced by `.abTEMFloat32`), and within 1e−6 with the kit's default exact π. Each slice has a written registration with its bounds and amendments:
  `docs/registration-slice-1.md` (potentials, multislice), `docs/registration-slice-2.md` (frozen phonons, scans,
  detectors), `docs/registration-slice-3.md` (CTF), `docs/registration-performance.md` (the Metal path) and
  `docs/registration-generality.md` (the defaults as parameters).
- **Truth tests** independent of abTEM: Kirkland's table, the k = 0 identity, the far field, vacuum propagation, slice
  symmetry, the probe, the aberration phase, the thickness series.
- **Mutants**: every new test is broken on purpose and must fail (`docs/slice-*-mutants-2026-10-15.md`).
- **Refuters**: an independent attempt to break the agreement on new inputs for each slice (recorded in the registrations).
- **`docs/parity.md`** lists abTEM's public API entry by entry: done, partial or absent (planned or not planned).
- **`tools/consumer-check`**: a separate package that uses the whole public API from outside, with negative controls.
- **`docs/status.md`** has the measured table of every gate.

## What is not claimed

Details in `docs/status.md`, "Known limits".

- Nothing is validated against a micrograph. Partial coherence is abTEM's envelope approximation; there is no detector
  response or noise. Frozen phonons use Einstein σ given by the caller, and agree with abTEM exactly only when abTEM's
  displacements are supplied (otherwise in distribution).
- Metal's float32 FFT is not exactly unitary (it loses about 8e−7 of the total intensity per slice; a per-step
  renormalisation brings that to 5e−7 or less); the Metal bounds are registered to 100 slices, and past 300 slices the
  pattern error exceeds its bound.
- Cells are orthogonal and always periodic; atoms in the gap of a short slice list are kept (abTEM drops them);
  `.afterSlices` always adds the last slice. Finite projections, rotated or cut cells and some of abTEM's options are not
  ported (`docs/parity.md`). `AtomicStructure` and the example's CIF reader are minimal local stand-ins.

## Running the example

```sh
swift run -c release SimulatoExample [file.cif] [output-directory]
```

Without arguments it uses the bundled `Examples/SimulatoExample/Al.cif`, builds a 3 × 3 × 10-cell crystal on 256², and writes
`Al-hrtem.png` (CTF at Scherzer defocus, 300 keV, Cs = 1 mm, focal spread 30 Å), `Al-haadf.png` (16 × 16 scan over one
cell, 70 to 200 mrad) and `Al-numbers.csv` to the current directory, printing the image statistics and the HAADF values,
on Metal when a device exists. It takes a few seconds.
The CIF must have an orthogonal cell.

## Tests

```sh
tools/test.sh    # swift test -c release; stamps .build/last-green on a clean tree
```

Never run `swift build --build-tests -c release`; `tools/test.sh` is the gate. The abTEM oracle (a Python venv) is not in the
repository: `docs/registration-slice-1.md` section 2 says how it is pinned, and `tools/slice-*/pins.py` regenerate the fixtures.

## Contributing and licence

Issues and pull requests are welcome. A change that moves a number needs a test that fails when the code is broken and a
comparison against abTEM; say which. GPL-3.0-or-later; see `LICENSE` and `NOTICE`.
