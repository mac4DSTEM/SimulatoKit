# Parity with abTEM

The coverage table of abTEM 1.0.10 (`References/python/abtem`, git `164e644`): its public classes and functions, their options, read from
the source in `site-packages/abtem` against SimulatoKit. Rewritten for 0.3.0 (the growing phase's kit audit); every row
has one state and one priority.

**Totals (128 rows): 37 present, 37 partial, 54 missing.** By priority: **study** 27 rows (11 present, 5 partial, 11 missing);
**useful** 79 (26, 32, 21); **not planned** 22 (all missing, each with its reason).

States:

- **present**: ported, the same numbers within the tolerance below, pinned by a test.
- **partial**: ported with options missing (named in the row), or ported and not yet pinned (named).
- **missing**: not ported; the row says "planned: <where>" or "not planned: <reason>".

Priority: **study** is needed for the TEM study the kit was built for and the apps that use it (`ROADMAP.md`), **useful**
is worth having, **not planned** has its reason in the row.

Tolerance. abTEM's wheel ships no tests (`tests/` is not in `site-packages`), so its own tolerances cannot be read from
disk. Where upstream states none the bound is `numpy.testing.assert_allclose`'s default for the
dtype: rtol 1e−7 for float64, 1e−6 for float32 (abTEM's complex64 centre of mass is held to 1e−5 of its largest value
because its float32 sum order is its own). The earlier slices pinned abTEM's exit waves and patterns at 1e−10 of the maximum
with `.abTEMFloat32`, and 1e−6 with the default `.exact` parameterisation; those pins are tighter than the default and are
kept. Readers, I/O, indices and metadata are exact.

Execution knobs (dask `lazy`, `max_batch`, `to_cpu`, `url`, `pbar`, `device`) are not counted as missing options; options that
change a result are.

Changes in 0.3.0 (rows marked 0.3.0): any grid length (was f · 2^m, f in 1, 3, 5, 15), `Potential(sampling:)`,
`Grid(…sampling…)`, `GridScan(…sampling…)`, beam tilt on `Waves`, `Probe` and the propagator (CPU and Metal), the centre-of-mass
detector, the Gaussian source size and Poisson noise on a scan.

## Beam, constants and grid

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `core.energy.energy2wavelength` | `ElectronBeam.wavelength` | **present** (bit-equal on the pins) | useful |
| `core.energy.energy2sigma` | `ElectronBeam.interactionParameter` | **present** (bit-equal on the pins) | useful |
| `core.energy.relativistic_mass_correction` | `ElectronBeam.relativisticMassFactor` | **partial**: ported, no test reads it (σ, pinned, carries the same factor; a pin of its own is missing) | useful |
| `core.energy.energy2mass` | | **missing**: planned: SimulatoCore (`ElectronBeam.mass` = the factor × mₑ) | useful |
| `Accelerator`, `HasAcceleratorMixin` (`energy`, `wavelength`, `sigma`, `lock_energy`, `match`, `check_match`) | `ElectronBeam(energy:)` | **partial**: energy, wavelength and σ; an undefined energy (None), `lock_energy` and `match` not offered; the energy must be finite and positive | useful |
| `core.constants.kappa` | `PhysicalConstants.kappa` | **present** (the same double) | useful |
| `Grid`, `HasGrid2DMixin` (`extent`, `gpts`, `sampling`, `reciprocal_space_sampling`, `endpoint`) | `Grid(nx:ny:extentX:extentY:)`, `Grid(extentX:extentY:samplingX:samplingY:)`, `.samplingX`/`Y`, `.reciprocalSamplingX`/`Y`, `.matches(_:)` | **partial**: any gpts ≥ 2 on both backends (0.3.0; vDSP where it has a length, Bluestein otherwise, pinned 90 x 75 against abTEM); `sampling=` gives ⌈extent / sampling⌉ points as abTEM; two dimensions, endpoint False only; the reciprocal sampling is abTEM's expression (bit-equal) | study |
| `core.grid.spatial_frequencies`, NumPy `fftfreq` | `Grid.fftfreq`, `.frequenciesX`/`Y` | **present** (exact today, `SimulatoCoreTests`) | useful |
| `BaseWaves.angular_sampling`, `.cutoff_angles`, `.antialias_cutoff_gpts`, `.antialias_valid_gpts`, `.rectangle_cutoff_angles`, `.full_cutoff_angles`, `.cutoff_frequencies` | `DiffractionIntensity.angularSampling`, `.cutoffAngles`, `Antialias.cutoffPatternSize(on:)` | **partial**: angular sampling, cutoff angles and the odd cutoff size; the valid size only inside `PixelatedDetector(.valid)`; the rectangle and full cutoff angles and the cutoff frequencies not offered | useful |
| `abtem.config` (dask, device, precision, antialias settings), `abtem.axes` (axis metadata) | | **missing**: not planned: global configuration and ensemble-axis metadata belong to abTEM's lazy array objects; the kit takes explicit arguments | not planned |
| `abtem.distributions` (parameter distributions for ensembles) | | **missing**: not planned: a caller loops over a parameter and averages (focal series, tilt series) | not planned |

## Potentials and slicing

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `Potential(atoms, gpts, slice_thickness, parametrization="lobato", projection="infinite")`, the lazy potential, and its `device` | `Potential(structure:nx:ny:sliceThickness:parameterization:)`, `.build() -> ProjectedPotential` | **present**: the grid, the slice thicknesses, the wrapped atoms with their slice and the parameterisation, without the arrays; `build()` equals the former eager construction (asserted exactly today); on `.metal` the potential is built on the GPU (V within 1e−5 × max V of the CPU's) | useful |
| `ase.Atoms` as the atoms argument (cell, positions, numbers, `pbc`) | `AtomicStructure(box:atoms:)`, `.repeated(_:_:_:)` | **partial**: orthogonal (diagonal) boxes, atomic number and Cartesian position only; periodic on all axes (abTEM honours ASE's `pbc`); repeats in ASE's order; no ASE or CIF reader in the library; a local type on purpose (the kit stays standalone, `docs/adr/0001-standalone-structure-types.md`) | useful |
| `Potential(slice_thickness=float or sequence)`, `SliceIndexedAtoms` | `sliceThickness:`, `sliceThicknesses:`, `Slicing.thicknesses`, `.validated`, `.wrapped`, `.sliceIndices` | **present**: ⌈Lz / s⌉ equal slices, or a list that must sum to Lz (NumPy `isclose`), periodic wrap, top-face rule; no default of 1 Å (the thickness is required); an atom in the gap of a short slice list is kept, abTEM drops it | useful |
| `Potential(sampling=)` | `Potential(structure:sampling:sliceThickness:)` (one value or an (x, y) pair) | **present**: ⌈box / sampling⌉ points per axis, pinned against abTEM's gpts on six samplings (0.3.0) | study |
| `Potential(exit_planes=)` | `ExitPlanes` on `Multislice.run` | **present**: `.every(m)` includes the entrance plane as abTEM does; `.afterSlices` always adds the last slice, abTEM's tuple does not | useful |
| `Potential(projection="finite")` (`QuadratureProjectionIntegrals`), `SlicedAtoms(xy_padding, z_padding)` | | **missing**: planned: SimulatoMultislice (finite projections for thick slices) | study |
| `Potential(periodic=False)` (cut cell with margins) | | **missing**: planned: SimulatoMultislice | study |
| `Potential(plane, origin, box)` (rotated and cut cells, `orthogonalize_cell`) | | **missing**: planned: with non-orthogonal cells in `AtomicStructure` (`plane="xy"`, `origin` 0 and the cell's own box only) | study |
| `Potential(integrator=)`, `ScatteringFactorProjectionIntegrals`, `GaussianProjectionIntegrals`, `ProjectionIntegralTable` | | **missing**: not planned: the infinite scattering-factor projection is the only integrator | not planned |
| `Potential.build(first_slice, last_slice, max_batch, lazy)`, `.generate_slices`, iteration | `Potential.build()` | **partial**: the whole potential only; `first_slice`, `last_slice`, `max_batch`, `lazy` (dask) and slice generation not offered | useful |
| `PotentialArray(array, slice_thickness, extent, sampling, exit_planes)` (the eager potential) | `ProjectedPotential` (`grid`, `slices`, `sliceOfAtom`, `init(grid:slices:)`) | **partial**: built by `Potential.build()` (done); an array from outside cannot be supplied: `PotentialSlice` has no public initialiser; no ensemble axes or stored exit planes | study |
| `PotentialArray.num_slices`, `.slice_thickness`, `.thickness`, `.slice_limits`, `.extent`, `.sampling` | `ProjectedPotential.slices.count`, `PotentialSlice.thickness`, `.thickness`, `.grid` | **partial**: `slice_limits` not offered | useful |
| `PotentialArray.transmission_function(energy)` | `Multislice.transmissionFunctions(_:beam:)` | **present**: the bandlimited t = IFFT(FFT(exp(iσV)) · A) | useful |
| `TransmissionFunction.transmit(waves, conjugate)` | | **missing**: planned: SimulatoMultislice (the transmission step on its own; it runs inside `Multislice.run`) | useful |
| `PotentialArray.project()`, `.to_images()` | | **missing**: planned: SimulatoMultislice (the summed slices as a `RealField`) | study |
| `PotentialArray.tile`, `CrystalPotential` (tiled potentials) | | **missing**: planned: SimulatoMultislice | study |
| `PotentialArray.depth_profile`, `.show`, `.show_depth_profile` | | **missing**: not planned: plotting (the kits have no UI) | not planned |
| `crystal_slice_thicknesses(atoms, tolerance)` | | **missing**: planned: SimulatoMultislice (`Slicing`) | study |
| `ChargeDensityPotential`, `GPAWPotential` | | **missing**: not planned: needs DFT inputs; outside the common use of the kit | not planned |
| `orthogonalize_cell`, `standardize_cell`, the helpers of `atoms.py` (`rotate_atoms`, `cut_cell`, `pad_atoms`, `shrink_cell`, ...) | | **missing**: planned: with DiffractoKit's cells | study |

## Parametrizations

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `LobatoParametrization` (`projected_scattering_factor`, `scaled_parameters`) | `LobatoParameterization` (`projectedScatteringFactor(atomicNumber:k2:)`, static `scaledParameters(atomicNumber:)`, `parameters(atomicNumber:)`, `supports(atomicNumber:)`) | **partial**: neutral atoms (`charge` not offered), Z = 1 to 103, the shipped table only (`parameters=` custom JSON and `sigmas=` not offered); `.abTEMFloat32` bit-close, `.exact` (the default) registered at 1e−7 | useful |
| `LobatoParametrization.potential`, `.projected_potential`, `.scattering_factor`, `.charge`, `.x_ray_scattering_factor`, `.finite_projected_potential`, `.line_profiles` | | **missing**: planned: SimulatoMultislice (the closed forms; `projected_potential` is the far-field truth's reference in the tests, evaluated in Python) | useful |
| `KirklandParametrization`, `PengParametrization`, `EwaldParametrization` | | **missing**: planned: SimulatoMultislice | study |
| `Parametrization.sigmas`, `sigma_to_B`, `B_to_sigma`, `validate_sigmas` (Debye–Waller tables) | | **missing**: planned: SimulatoMultislice (σ is given by the caller per element today) | study |
| `Parametrization(parameters)`, `.to_json`, `.from_json`, `.fit` | | **missing**: not planned: custom tables and fitting are research tools | not planned |

## Waves and multislice

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `PlaneWave(energy, extent, gpts, sampling, normalize=False, device)`, `.build()` | `Waves.planeWave(beam:grid:tilt:)` | **partial**: `normalize=True` not offered; the grid is given as gpts and extent (or built from a sampling by `Grid`); the backend is chosen at `Multislice.run` | useful |
| `PlaneWave(tilt=)`, `Probe(tilt=)`, `Waves.base_tilt` | `Waves.planeWave(beam:grid:tilt:)`, `Probe(…tilt:)`, `Waves.tilt` [mrad] | **present**: as abTEM the tilt is metadata that enters the Fresnel propagator, e^{−2πi (kx tan tx + ky tan ty) Δz}; exit waves keep it; pinned on a 90 x 75 grid (plane wave and probe), CPU and Metal; a `.metal` batch must share one tilt (0.3.0). Tilt ensembles (distributions of tilts) not offered | study |
| `Waves(array, energy, extent, sampling, reciprocal_space, ensemble_axes_metadata, metadata)` | `Waves(field:grid:beam:normalization:tilt:)` | **partial**: one two-dimensional wave in real space with its tilt; no ensemble axes, reciprocal-space storage or other metadata; `normalization` is `.values` or `.reciprocalSpace` | useful |
| `Probe(semiangle_cutoff, energy, extent, gpts, sampling, soft=True, aperture, aberrations)` | `Probe(beam:grid:aperture:aberrations:tilt:)` | **partial**: `tilt` offered (0.3.0); no `scan_positions`, `metadata`, distributions of parameters, or aberration coefficients as keywords (set on `Aberrations`); a probe is rebuilt to change a parameter | useful |
| `Probe.build(scan=)` | `Probe.build(at:)` | **partial**: one position per call (a scan: `STEM.scan`); `scan=` as a grid of positions not offered | useful |
| `Probe.scan(potential, scan, detectors)`, `Probe.multislice(potential, scan, detectors)` | `STEM.scan(probe:potentials:scan:detectors:backend:)` | **present**: positions bit-equal, detectors averaged over the potentials; `max_batch` is chosen inside (chunks of at most 256 MB of waves), `lazy` not offered | useful |
| `Waves.multislice(potential, detectors)`, `PlaneWave.multislice(...)` | `Multislice.run(_:through:exitPlanes:backend:)` | **partial**: waves in, waves at the exit planes out; no `detectors=` (apply a `Detector` to `DiffractionIntensity(exit)`, as `STEM.scan` does per position); a `WavesDetector` over a scan is `STEM.scan(... detector:)` | useful |
| `fft_interpolate(array, new_shape)` for waves | `Waves.resampled(rows:columns:)` | **present**: crop or pad in Fourier space, `normalization="values"`; the 2D case only; pinned through the WavesDetector rows | study |
| `multislice_and_detect(waves, potential, algorithm=FourierMultislice(order=1))` | `Multislice.run(_:through:exitPlanes:backend:)` | **present**: `.cpu` float64 and `.metal` float32 (Backends) | useful |
| `FourierMultislice(order=2)`, `conjugate`, `transpose` | | **missing**: planned: SimulatoMultislice | useful |
| `RealSpaceMultislice` (finite differences, backscattering) | | **missing**: not planned for now (a research option) | not planned |
| `FresnelPropagator.get_array(waves, thickness, order=1)` | `Multislice.propagator(grid:beam:thickness:tilt:)` | **present**: order 1, bandlimited, with the beam tilt (0.3.0); pinned through the exit waves | useful |
| `FresnelPropagator.propagate(waves, thickness, in_place, order)` | | **missing**: planned: SimulatoMultislice (free propagation of a wave over Δz; order 2 with `FourierMultislice(order=2)`) | useful |
| `AntialiasAperture`, `antialias_aperture(gpts, sampling)`, `bandlimit` | `Antialias.aperture(on:)`, `.cutoff`, `.taper` | **present**: abTEM's defaults 0.6666666 and 0.01 (not configurable) | useful |
| `Waves.intensity()` | `Waves.intensity()` | **present** | useful |
| `Waves.diffraction_patterns(max_angle="cutoff" or "full", fftshift=True)` | `Waves.diffractionPattern(maxAngle:)` | **partial**: `"valid"`, a numeric angle, `block_direct`, `return_complex`, `renormalize=False`, `fftshift=False` and `parity` other than "odd" not offered | useful |
| `Waves.phase()`, `.real()`, `.imag()`, `.to_images()` | `ComplexField.real`, `.imaginary`, `.squaredMagnitude()` | **partial**: raw arrays; no phase helper and no image-returning methods | useful |
| `Waves.apply_ctf(ctf)` (HRTEM image formation) | `Waves.applying(_:)` | **present**: the `**kwargs` shorthand that builds a CTF is not offered | useful |
| `Waves.normalize`, `.tile`, `.downsample`, `.convolve`, `.phase_shift`, `.ensure_reciprocal_space`, `.ensure_real_space` | | **missing**: planned as needed: SimulatoCore | useful |
| `Waves.scan`, `.depth_profile`, `.show_depth_profile`, `.show`, `Probe.profiles`, `.ctf`, `.show` | | **missing**: not planned: plotting and lazy-array helpers | not planned |
| `Waves.transition_potential_multislice`, `Probe.transition_potential_scan` | | **missing**: not planned: inelastic scattering (see Other abTEM capabilities) | not planned |

## Probe and transfer

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `Aperture(semiangle_cutoff, soft=True)`, `soft_aperture`, `hard_aperture` | `Aperture(semiangleCutoff:soft:)`, `.values(on:wavelength:)` | **present** | useful |
| `Aberrations` coefficients C10 to C56 and φ (25 symbols), χ(α, φ) | `Aberrations` (`C10` ... `phi56`), `.chi(alpha:phi:wavelength:)` | **present**: orders 1 to 5 | useful |
| `Aberrations` aliases: `defocus` = −C10, `Cs` = C30, and 23 more (`C5`, `astigmatism`, `astigmatism_angle`, `coma`, `trefoil`, `quadrafoil`, `pentafoil`, `hexafoil`, `astigmatism3`, `astigmatism5`, `coma4`, `trefoil4`, `quadrafoil5` and their `_angle`) | `Aberrations.defocus`, `.Cs` | **partial**: `defocus` and `Cs` only; the other 23 by symbol | useful |
| `Aberrations(aberration_coefficients=dict, **kwargs)`, `set_aberrations`, `.aberration_coefficients` | | **missing**: not planned: the coefficients are stored properties, set directly | not planned |
| `CTF(semiangle_cutoff, soft, focal_spread, angular_spread, aberration_coefficients, energy, extent, gpts, sampling, flip_phase, wiener_snr)` | `ContrastTransferFunction(aberrations:semiangleCutoff:softAperture:focalSpread:angularSpread:)`, `.kernel(on:wavelength:)` | **partial**: `flip_phase` and `wiener_snr` not offered; energy and grid are given at `kernel(on:wavelength:)`, not stored; negative spreads refused (abTEM overflows) | useful |
| `SpatialEnvelope(angular_spread, aberration_coefficients)`, `TemporalEnvelope(focal_spread)` (partial coherence) | `ContrastTransferFunction.spatialEnvelope`, `.temporalEnvelope` | **present**: values at (α, φ) and inside the kernel; the quasi-coherent approximation, as abTEM | useful |
| `transfer.scherzer_defocus`, `point_resolution`, `nyquist_sampling` (functions) | `ContrastTransferFunction.scherzerDefocus`, `.pointResolution`, `.nyquistSampling` | **present** | useful |
| `CTF.scherzer_defocus`, `.point_resolution`, `.crossover_angle` (from the CTF's own Cs and energy), `BaseAperture.nyquist_sampling` | | **missing**: planned: SimulatoProbe (instance forms; the statics take Cs and energy) | useful |
| `CTF.to_point_spread_functions`, `.profiles`, `.to_diffraction_patterns` | | **missing**: planned: SimulatoProbe | useful |
| `CTF.show`, `BaseTransferFunction.show` | | **missing**: not planned: plotting | not planned |
| `Bullseye`, `Vortex`, `AnnularAperture`, `Zernike`, `RadialPhasePlate` | | **missing**: planned as needed: SimulatoProbe | useful |
| `transfer.polar2cartesian`, `cartesian2polar` | | **missing**: planned: SimulatoProbe | useful |

## Frozen phonons

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `FrozenPhonons(atoms, num_configs, sigmas, directions, ensemble_mean, seed)` | `FrozenPhonons(structure:sigmas:configurations:seed:directions:)`, `FrozenPhonons(structure:displacements:)` | **partial**: see the five rows below | useful |
| `sigmas` | `sigmas: [Int: Double]` (Å, per atomic number) | **partial**: a σ for every element of the structure, as abTEM requires; a scalar for all atoms, symbols as keys, per-atom arrays and anisotropic σ (a 3-vector) not offered | useful |
| `num_configs` | `configurations:`, `FrozenPhonons.count` | **present** | useful |
| `seed` | `seed: UInt64`, `FrozenPhonons.configurationSeed(_:_:)` | **partial**: abTEM's displacements reproduced only when given (`displacements:`, bit-equal); NumPy's PCG64 normals are not ported, the kit's generator (`GaussianGenerator`: SplitMix64, xoshiro256**, Box–Muller) agrees in distribution (KS tests); a seed is required (no None) | useful |
| `directions="xyz"` | `FrozenPhonons.Directions` (`.x`, `.y`, `.z`, `.all`) | **present** | useful |
| `ensemble_mean=True` | averaging inside `STEM.scan` | **partial**: the mean only; `ensemble_mean=False` (one measurement per configuration) not offered: planned: SimulatoDetector (a caller can scan each potential alone) | useful |
| `FrozenPhonons.atoms`, `.randomize`, `.num_configs`, `Potential(frozen_phonons)` | `FrozenPhonons.configuration(_:)`, `.displacements`, `.potentials(nx:ny:sliceThickness:parameterization:)` (lazy: `lazyPotentials`) | **present**: one potential per configuration, averaged by `STEM.scan` | useful |
| `FrozenPhonons.to_atoms_ensemble`, `AtomsEnsemble(trajectory)`, `DummyFrozenPhonons` | | **missing**: planned: SimulatoMultislice (a trajectory of structures, for example from molecular dynamics) | useful |

## Scans

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `GridScan(start, end, gpts, endpoint=False)` | `GridScan(start:end:nx:ny:endpoint:)`, `.positions`, `GridScan.linspace` | **present**: positions bit-equal to abTEM's (NumPy's `linspace`, x the first index); a scan of extent ≤ 0 on both axes refused, as abTEM | useful |
| `GridScan(endpoint=True)` | `endpoint: true` | **partial**: implemented as `linspace(endpoint: true)`; the pins use `endpoint=False`, so it is not pinned against abTEM | study |
| `GridScan(sampling=)` | `GridScan(start:end:samplingX:samplingY:endpoint:)` | **present**: ⌈extent / sampling⌉ points, pinned on three scans against abTEM (0.3.0); abTEM's default sampling from the probe's Nyquist (`_validate_scan_sampling`) not offered | study |
| `GridScan(potential=, fractional=, end=None)` (extent from the potential, fractional coordinates, `Atom` start and end) | | **missing**: planned: SimulatoProbe (`start` and `end` are required, in Å) | study |
| `LineScan(start, end, gpts, sampling, endpoint=True, fractional, potential)`, `.at_position`, `.add_margin` | `LineScan(start:end:gpts:sampling:endpoint:fractional:extent:)`, `.addingMargin`, `LineScan.atPosition`, `.positions`, `.shape` | **present**: pinned on 13 lines against abTEM (positions bit-equal, the end points to the last bit after a margin); a sampling wins over gpts and a margin recomputes the count end by end, as abTEM's setters do; `end` and one of `gpts`/`sampling` are required (abTEM's `end=None` and `match_probe` defaults not offered); `potential=` is the extent `extent:`; `STEM.scan` over a line: 5e-15 of the ADF with abTEM's float32 π, 1.1e-7 with the exact π, `.metal` 1.5e-6 of `.cpu` | study |
| `CustomScan(positions)` | | **missing**: planned: SimulatoProbe (`Probe.build(at:)` takes any position today; `STEM.scan` takes a `GridScan`) | study |
| `BaseScan.get_positions`, `.limits`, `.num_positions`, `.shape` | `GridScan.positions`, `.nx`, `.ny`, `LineScan.positions`, `.gpts`, `ScanPattern.shape` | **partial**: positions and the shape; `limits` not offered | useful |
| `GridScan.add_to_plot`, `LineScan.add_to_plot` | | **missing**: not planned: plotting | not planned |

## Detectors

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `AnnularDetector(inner, outer, offset)` | `AnnularDetector(inner:outer:)` | **partial**: `offset` not offered; a negative inner is refused (abTEM squares it); `to_cpu`, `url` do not apply | useful |
| `PixelatedDetector(max_angle, resample, reciprocal_space)` | `PixelatedDetector(maxAngle:)` | **partial**: "cutoff", "valid", "full"; numeric angles, `resample` and real-space detection not offered | useful |
| `SegmentedDetector(nbins_radial, nbins_azimuthal, inner, outer, rotation, offset)` | `SegmentedDetector(radialBins:azimuthalBins:inner:outer:rotation:)` | **partial**: summed in float64 where abTEM sums in float32; inner = outer refused; `offset` not offered | useful |
| `FlexibleAnnularDetector(step_size, inner, outer)` | | **missing**: planned: SimulatoDetector | useful |
| `WavesDetector(gpts)` | `WavesDetector(rows:columns:)`, `STEM.scan(probe:potentials:scan:detector:)` returning `ScanWaves` | **present**: the exit wave per position, `gpts` as a Fourier crop or pad; pinned at 1.8e-15 of the maximum on a 3 x 2 scan, six `gpts` (crop, pad, mixed, odd, even), `.metal` at 3.6e-7 against `.cpu`; one potential only (a frozen-phonon list is refused: averaging complex waves is not abTEM's `ensemble_mean` of intensities; scan each configuration alone); refuses above half of the physical memory; `to_cpu`, `url` do not apply | study |
| `BaseDetector.detect(waves)`, `.angular_limits` | `Detector.detect(_:)` on a `DiffractionIntensity`, `.shape(for:beam:)` | **present** | useful |
| `AnnularDetector.get_detector_region`, `get_detector_regions`, `.show` | | **missing**: not planned: plotting | not planned |

## Measurements

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `Images(array, sampling)` | `RealField` (`nx`, `ny`, `values`, `samplingX`, `samplingY`, `sum`, `max`), `ScanMeasurement.image(samplingX:samplingY:)` | **partial**: a two-dimensional array with its sampling; no ensemble axes, metadata or units | useful |
| `DiffractionPatterns(array, sampling, fftshift)` | `RealField` from `Waves.diffractionPattern`, `DiffractionIntensity` (unshifted, with `angularSampling`) | **partial**: the pattern with its sampling in 1/Å; no angular sampling or units stored on the `RealField` | useful |
| `PolarMeasurements` (`integrate`, `integrate_radial`, `differentials`, `to_diffraction_patterns`) | `ScanMeasurement.integrate(detectorRegions:)`, `.integrate(detector:)`, `.integrateRadial(inner:outer:detector:)`, `.differentials(direction1:direction2:)` on a `SegmentedDetector` scan | **partial**: regions, radial limits (truncated to bins as abTEM) and the DPC differentials (second regions minus first, abTEM's order) against abTEM on a planted scan, float64 sums at rtol 1e-7, differentials at 1e-6 (abTEM complex64); the differentials come as two scans, not one complex image; `azimuthal_limits` not offered (abTEM divides them by the radial sampling), `to_diffraction_patterns` not offered; refused where abTEM is silent: no region (abTEM: zeros), a region or radial index below 0 (Python reads it from the end) | useful |
| `RealSpaceLineProfiles`, `ReciprocalSpaceLineProfiles`, `_BaseMeasurement1D.width`, `.interpolate` | `RealSpaceLineProfiles`, `ReciprocalSpaceLineProfiles`, `LineProfiles.width(height:)`, `.interpolate(sampling:gpts:order:endpoint:)`, `ScanMeasurement.lineProfiles(for:)` | **present**: width to 4.5e-8 and interpolate to 5.5e-8 of the maximum against abTEM (abTEM stores both as float32; `assert_allclose` rtol 1e-7 holds at 0.6 of its bound), on the scan's profile and on planted profiles with ensemble axes, orders 0 to 5, both endpoints; `lineProfiles` only for detectors of one number or a vector (abTEM returns other measurement types for segmented and pixelated); `RealSpaceLineProfiles.tile`, `ReciprocalSpaceLineProfiles.angular_extent`, `show` and the lazy path not offered | useful |
| `MeasurementsEnsemble`, `BaseMeasurements.to_measurement_ensemble`, `normalize_ensemble` | | **missing**: not planned: Swift arrays; a caller keeps the configurations apart | not planned |
| `IndexedDiffractionPatterns`, `DiffractionPatterns.index_diffraction_spots` | | **missing**: planned: with DiffractoKit's indexing | useful |
| `BaseMeasurements.intensity()` | `Waves.intensity()`, `ComplexField.squaredMagnitude()` | **present** | useful |
| `BaseMeasurements.real()`, `.imag()` | `ComplexField.real`, `.imaginary` | **partial**: raw arrays, not measurements | useful |
| `BaseMeasurements.phase()`, `.abs()` | | **missing**: planned: SimulatoCore | useful |
| `BaseMeasurements.reduce_ensemble()` | the mean in `STEM.scan` | **partial**: the mean over configurations only | useful |
| `BaseMeasurements.relative_difference(other, min_relative_tol)` | | **missing**: planned: SimulatoCore | useful |
| `poisson_noise(dose_per_area, total_dose, samples, seed)` (`BaseMeasurements`, `DiffractionPatterns`, `PolarMeasurements`) | `ScanMeasurement.poissonNoise(dose:seed:)`, `.poissonNoise(dosePerArea:samplingX:samplingY:seed:)` | **partial**: values × dose, clipped at 0, one Poisson draw each; mean and variance pinned statistically from 0.5 to 800 electrons; the kit's own generator, not NumPy's stream (a seed repeats on the kit, not against abTEM); `samples > 1` (an ensemble axis) and noise on an `Images` / `Waves` object not offered (0.3.0) | study |
| `DiffractionPatterns.integrate_radial(inner, outer, offset)` | `AnnularDetector.detect(_:)` | **partial**: only as a detector on an exit wave's `DiffractionIntensity`; not a method of a stored pattern; `offset` not offered | useful |
| `DiffractionPatterns.polar_binning(nbins_radial, nbins_azimuthal, inner, outer, rotation, offset)` | `SegmentedDetector.detect(_:)` | **partial**: as the row above | useful |
| `DiffractionPatterns.radial_binning(step_size, inner, outer)`, `.azimuthal_average` | | **missing**: planned: SimulatoDetector | useful |
| `DiffractionPatterns.center_of_mass(units)` | `CenterOfMassDetector(maxAngle:units:)` | **partial**: Σ P·α over the pixelated detector's crop in mrad or 1/Å, not divided by Σ P (as abTEM), pinned within 1e−5 of abTEM's complex64 result (0.3.0); a detector on the exit wave, not a method of a stored pattern; `integrated_center_of_mass` is `ScanMeasurement.integratedCenterOfMass` (next row) | study |
| `DiffractionPatterns.integrated_center_of_mass()` | `ScanMeasurement.integratedCenterOfMass(samplingX:samplingY:)` of a `CenterOfMassDetector` scan | **present**: abTEM's centre of mass through the kit 2e−16 (rtol 1e−7 default); the kit's own centre of mass 5e−8 (its COM is registered at 1e−5); no scaling beyond abTEM's (units of the COM times Å, minimum 0) | study |
| `Images.integrate_gradient` | `RealField.integratingGradient(x:y:)`, `ComplexField.integratedGradient(samplingX:samplingY:)` | **present** within assert_allclose float64 defaults (rtol 1e−7, atol 1e−12 of the max); measured 8e−16 of the maximum on even, odd, mixed, 2×3 and random (Nyquist) fields; planted potentials recovered to 1e−9 of their range | study |
| `integrate_disc` (`Images.integrate_disc`) | `RealField.integrateDisc(x:y:radius:mean:border:)` | **present** same tolerance; measured 2e−15 relative; `interpolate` not offered | study |
| `Images.tile` | `RealField.tiled(repetitionsX:repetitionsY:)` | **present** exact | useful |
| `DiffractionPatterns.crop(max_angle, max_frequency, gpts)` | `Waves.diffractionPattern(maxAngle: .cutoff / .full)`, `PixelatedDetector(maxAngle:)` | **partial**: the crop is chosen when the pattern is made; not a method of a stored pattern | useful |
| `DiffractionPatterns.interpolate(sampling, gpts)`, `Images.interpolate(sampling, gpts, method, boundary, order, normalization, cval)` | | **missing**: planned: SimulatoCore | useful |
| `DiffractionPatterns.gaussian_source_size(sigma)` | `ScanMeasurement.gaussianSourceSize(sigma:samplingX:samplingY:)` | **present**: scipy `gaussian_filter` over the two scan axes, mode wrap, radius ⌊4σ + ½⌋; pinned at rtol 1e−7 (float64) incl. a σ that wraps the scan (0.3.0) | study |
| `DiffractionPatterns.bandlimit`, `.block_direct`, `.tile_scan` | | **missing**: planned: SimulatoDetector | useful |
| `Images.crop(extent, offset, centered)` | | **missing**: planned: SimulatoCore | useful |
| `_BaseMeasurement2D.gaussian_filter(sigma, boundary, cval)` | | **missing**: planned: SimulatoCore | useful |
| `interpolate_line`, `interpolate_line_at_position` | `RealField.interpolateLine`, `.interpolateLineAtPosition` (images), `.interpolateReciprocalLine`, `.interpolateReciprocalLineAtPosition` (centred patterns) | **present**: 38 lines on a 24 × 18 image and a 21 × 20 pattern against abTEM, orders 0 to 5, width, margin, fractional, both endpoints: 3.5e-14 of the maximum, `assert_allclose` ratio below 1.2e-5 of its bound; `end` is required (abTEM's `(0, extent[0])` default not offered); a pattern is assumed centred (zero frequency at n / 2); complex and lazy arrays not offered; in SimulatoDetector, not SimulatoCore, because `LineScan` lives in SimulatoProbe | study |
| `Images.scan_noise(dwell_time, flyback_time, rms_power, ...)`, `Images.diffractograms()` | | **missing**: not planned for now: scan-noise modelling; diffractograms are the field of an image-analysis kit | not planned |
| `BaseMeasurements.show`, `to_zarr`, `to_tiff`, `to_hyperspy`, `to_data_array`, `to_quantem` | | **missing**: not planned: plotting and file formats | not planned |

## Other abTEM capabilities

| abTEM | SimulatoKit | State | Priority |
|---|---|---|---|
| `SMatrix`, `SMatrixArray` (PRISM) | | **missing**: not planned in 2026 (a later acceleration step) | not planned |
| `BlochWaves`, `StructureFactor` | | **missing**: not planned for now (DiffractoKit owns structure factors; Bloch waves a possible later slice) | not planned |
| `inelastic.core_loss`, `inelastic.plasmons`, `transition_potential_*` | | **missing**: not planned: SpectroKit's domain if ever | not planned |
| `mtf.py` (detector MTF), `noise.py` | | **missing**: planned: SimulatoDetector (with the noise rows above) | useful |
| `reconstruct.py` (ptychography), `mcf.py`, `magnetism`, `finite_difference.py` | | **missing**: not planned: research modules outside the kit's scope | not planned |
| `array.stack`, `concatenate`, `from_zarr`, dask laziness | | **missing**: not planned: Swift arrays; laziness is FormatoKit's | not planned |
| `visualize`, `show_atoms` | | **missing**: not planned: the kits have no UI | not planned |

## Backends

abTEM runs on NumPy, or on CuPy with `device="gpu"`. SimulatoKit has two backends behind one API (`Backend`, passed to
`Multislice.run`, `Multislice.runBatch` and `STEM.scan`; a `Potential` on `.metal` builds on the GPU):

| | Precision and engine | Role | Agreement |
|---|---|---|---|
| `.cpu` | float64, Accelerate; a batch of positions on at most 8 threads, equal to sequential runs (exact today) | the reference | against abTEM 1.0.10 (float64 configuration) at the pins' bounds: the potential and multislice pins (V 1e−11 × max V, exit waves, intensities and patterns 1e−10 × max with `.abTEMFloat32`; 1e−6 with `.exact`), the detector pins (annular and pixelated 1e−10, segmented 1e−6) and the CTF pins (kernels 1e−12, images 1e−10) |
| `.metal` | float32, MPSGraph FFT | abTEM's `device="gpu"` counterpart | against `.cpu`, `MetalAccuracy.registered`: exit wave 1e−4, intensity and pattern 2e−4 of the maximum; Σ\|ψ\|² 1.5e−6 per slice; scans 2e−4; V built on the GPU 1e−5 × max V, GPU probe 1e−5 × max \|ψ\|, grids to 1024²; a 90 × 75 grid with and without tilt (0.3.0) |

The Metal bounds hold to 100 slices; MPSGraph's float32 FFT is not exactly unitary, a known
limit (`tools/metal-fft-gain`; README "Known limits"). The CTF (`Waves.applying`) runs on the CPU only. A Metal result beyond a bound is a
failure of the Metal path, never a tolerance to widen. `MetalAccuracy.registered` carries these bounds as a value.

Execution choices are parameters (`ExecutionOptions`, passed as `options:`; 0.2.0):
threads, Metal batch memory and shapes, the CPU scan chunk, the renormalisation, the executable cache. abTEM's nearest
counterparts are its `max_batch` and dask settings, which are not counted as missing options above.

## Not abTEM's, SimulatoKit's own

- The Metal backend (`Backend.metal`) and the batched scan (`Multislice.runBatch`), registered against the CPU path; MPSGraph's float32 FFT gain is a known limit
  (`tools/metal-fft-gain`).
- The CPU batch on up to `ExecutionOptions.maximumCPUThreads` (8 by default) threads, equal to sequential runs (exact today).
- The constants of the CTF closed forms as parameters (`factor:` of `scherzerDefocus`, 1.5 by default as abTEM; `divisor:` of
  `pointResolution`, 6; `oversampling:` of `nyquistSampling`, 4).
- `LobatoParameterization.PiConstant`: the exact π by default, abTEM's float32 π on request.