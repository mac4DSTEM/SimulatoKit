# Parity with abTEM

The checklist of abTEM 1.0.10's documented public API (`abtem/__init__.py` at `164e644`, its `__all__` and the classes
behind it) against SimulatoKit, so "complete" is a measured claim (suite `ROADMAP.md`, "Kits", phase 4). Updated at every
tag and at every slice. States: **ported** (same numbers within the registered tolerance, pinned by tests), **partial**
(the default path ported, options not), **planned** (a named slice), **not planned** (with the reason).

State at 2026-10-15, slice 1 (`docs/registration-slice-1.md`): 17 entries ported or partial, 19 planned, 7 not planned.

## Beam and constants

| abTEM | SimulatoKit | State |
|---|---|---|
| `core.energy.energy2wavelength` | `ElectronBeam.wavelength` | ported (bit-equal on the pins) |
| `core.energy.energy2sigma` | `ElectronBeam.interactionParameter` | ported (bit-equal on the pins) |
| `core.energy.relativistic_mass_correction` | `ElectronBeam.relativisticMassFactor` | ported |
| `core.constants.kappa` | `PhysicalConstants.kappa` | ported (the same double) |

## Potentials

| abTEM | SimulatoKit | State |
|---|---|---|
| `Potential(atoms, gpts, sampling, slice_thickness, parametrization="lobato", projection="infinite")` | `ProjectedPotential(structure:nx:ny:sliceThickness:parameterization:)` | partial: `gpts` only (no `sampling` argument), orthogonal cells only, `plane="xy"`, `origin` 0, no `box`; periodic on all axes (abTEM honours ASE's `pbc`: amendment 3, item 5); an atom in the gap of a short slice list is kept, abTEM drops it (item 3) |
| `Potential(..., projection="finite")` (`QuadratureProjectionIntegrals`) | | planned (a later slice: finite projections for thick slices) |
| `Potential(..., exit_planes=)` | `ExitPlanes` on `Multislice.run` | ported (`.every(m)` includes the entrance plane as abTEM does; `.afterSlices` always adds the last slice, abTEM's tuple does not: amendment 3, item 6) |
| `Potential(..., periodic=False)` (cut cell with margins) | | planned |
| `Potential(..., plane, origin, box)` (rotated and cut cells, `orthogonalize_cell`) | | planned (after DiffractoKit's cell types) |
| `LobatoParametrization` (`projected_scattering_factor`, `scaled_parameters`) | `LobatoParameterization` | ported (`.abTEMFloat32` bit-close; `.exact` default, registered 1e−7) |
| `LobatoParametrization.potential`, `.projected_potential`, `.scattering_factor`, `.charge`, `.x_ray_scattering_factor` | | planned (the closed forms; `projected_potential` is the far-field truth's reference in the tests, evaluated in Python) |
| `KirklandParametrization`, `PengParametrization`, `EwaldParametrization`, `sigmas` (Debye–Waller) | | planned (slice 2 with frozen phonons for the sigmas) |
| `PotentialArray`, `.transmission_function(energy)` | `ProjectedPotential.slices`, `Multislice.transmissionFunctions` | ported |
| `CrystalPotential` (tiled potentials) | | planned |
| `ChargeDensityPotential`, `GPAWPotential` | | not planned (needs DFT inputs; outside the field's common use of the suite) |
| `atoms.orthogonalize_cell`, `standardize_cell` | | planned (with DiffractoKit's cells) |

## Waves and multislice

| abTEM | SimulatoKit | State |
|---|---|---|
| `PlaneWave(energy, normalize=False)` | `Waves.planeWave(beam:grid:)` | ported (`normalize=True` not offered) |
| `PlaneWave(tilt=)`, `Probe(tilt=)` | | planned |
| `Probe(semiangle_cutoff, energy, aberrations, soft=True)`, `.build(scan=position)` | `Probe(beam:grid:aperture:aberrations:)`, `.build(at:)` | ported (one position; distributions of parameters not offered) |
| `Probe.build(scan=GridScan/LineScan)`, `CustomScan`, `GridScan`, `LineScan` | | planned (slice 2, with the detectors) |
| `multislice_and_detect`, `FourierMultislice(order=1)` | `Multislice.run(_:through:exitPlanes:backend:)` | ported, `.cpu` float64 and `.metal` float32 |
| `FourierMultislice(order=2)`, `conjugate`, `transpose` | | planned |
| `RealSpaceMultislice` (finite differences, backscattering) | | not planned for now (a research option) |
| `Waves.intensity()` | `Waves.intensity()` | ported |
| `Waves.diffraction_patterns(max_angle="cutoff" or "full", fftshift=True)` | `Waves.diffractionPattern(maxAngle:)` | ported (`"valid"`, a numeric angle, `block_direct`, `return_complex`, `parity` other than "odd" not offered) |
| `Waves.phase()`, `.real()`, `.imag()`, `.complex_images()` | `ComplexField` (real, imaginary) | partial (raw arrays; no phase helper) |
| `Waves.apply_ctf(ctf)` (HRTEM image formation) | | planned (with `CTF`) |
| `Waves.normalize()`, `.tile()`, `.interpolate()` | | planned as needed |

## Transfer

| abTEM | SimulatoKit | State |
|---|---|---|
| `Aperture(semiangle_cutoff, soft=True)` | `Aperture(semiangleCutoff:soft:)` | ported |
| `Aberrations` (C10 ... C56, φ, the aliases `defocus` = −C10, `Cs`) | `Aberrations` | ported (aliases `defocus`, `Cs`; the others by symbol) |
| `CTF` | | planned (HRTEM images from a plane-wave exit) |
| `SpatialEnvelope`, `TemporalEnvelope` (partial coherence) | | planned |
| `Bullseye`, `Vortex`, `AnnularAperture`, `Zernike`, `RadialPhasePlate` | | planned as needed |
| `transfer.scherzer_defocus`, `point_resolution`, `nyquist_sampling`, `polar2cartesian` | | planned (small closed forms) |

## Detectors and measurements (slice 2)

| abTEM | SimulatoKit | State |
|---|---|---|
| `AnnularDetector`, `FlexibleAnnularDetector`, `SegmentedDetector`, `PixelatedDetector`, `WavesDetector` | `SimulatoDetector` | planned (slice 2) |
| `FrozenPhonons`, `AtomsEnsemble` | | planned (slice 2) |
| `Images`, `DiffractionPatterns`, `PolarMeasurements`, line profiles (`measurements.py`) | `RealField` | partial (arrays with their sampling; the measurement operations are planned with the detectors) |
| `IndexedDiffractionPatterns` | | planned (with DiffractoKit's indexing) |

## Other abTEM capabilities

| abTEM | SimulatoKit | State |
|---|---|---|
| `SMatrix`, `SMatrixArray` (PRISM) | | not planned in 2026 (a later acceleration lane) |
| `BlochWaves`, `StructureFactor` | | not planned for now (DiffractoKit owns structure factors; Bloch waves a possible later slice) |
| `inelastic.core_loss`, `inelastic.plasmons` | | not planned (SpectroKit's domain if ever) |
| `array.stack`, `concatenate`, `from_zarr`, dask laziness | | not planned (Swift arrays; laziness is FormatoKit's) |
| `visualize`, `show_atoms` | | not planned (the kits have no UI) |

## Not abTEM's, SimulatoKit's own

- The Metal backend (`Backend.metal`), registered against the CPU path; MPSGraph's float32 FFT gain is a known limit
  (registration amendment 2).
- `LobatoParameterization.PiConstant`: the exact π by default, abTEM's float32 π on request.
