# Parity with abTEM

The checklist of abTEM 1.0.10's documented public API (`abtem/__init__.py` at `164e644`: its `__all__`, the classes behind
it, their public methods and arguments, read from the docstrings: the pinned checkout carries no `docs/`) against
SimulatoKit, so "complete" is a measured claim (suite `ROADMAP.md`, "Kits", phase 4). Updated at every tag and at every
slice. Every row has exactly one state:

- **done**: ported, the same numbers within a registered tolerance, pinned by tests.
- **partial**: ported with options or arguments missing (named in the row), or ported and not yet pinned (named).
- **absent**: not ported; the row says "planned: <where>" or "not planned: <reason>".

Execution knobs (dask `lazy`, `max_batch`, `to_cpu`, `url`, `pbar`, `device`) are not counted as missing options; options that
change a result are.

Count at 2026-10-15, slices 1 to 3 and the performance lane (`docs/registration-slice-<n>.md`,
`docs/registration-performance.md`): **24 done, 34 partial, 64 absent** (122 rows).

## Beam, constants and grid

| abTEM | SimulatoKit | State |
|---|---|---|
| `core.energy.energy2wavelength` | `ElectronBeam.wavelength` | **done** (bit-equal on the pins) |
| `core.energy.energy2sigma` | `ElectronBeam.interactionParameter` | **done** (bit-equal on the pins) |
| `core.energy.relativistic_mass_correction` | `ElectronBeam.relativisticMassFactor` | **partial**: ported, no test reads it (σ, pinned, carries the same factor; a pin of its own is missing) |
| `core.energy.energy2mass` | | **absent**: planned: SimulatoCore (`ElectronBeam.mass` = the factor × mₑ) |
| `Accelerator`, `HasAcceleratorMixin` (`energy`, `wavelength`, `sigma`, `lock_energy`, `match`, `check_match`) | `ElectronBeam(energy:)` | **partial**: energy, wavelength and σ; an undefined energy (None), `lock_energy` and `match` not offered; the energy must be finite and positive |
| `core.constants.kappa` | `PhysicalConstants.kappa` | **done** (the same double) |
| `Grid`, `HasGrid2DMixin` (`extent`, `gpts`, `sampling`, `reciprocal_space_sampling`, `endpoint`) | `Grid(nx:ny:extentX:extentY:)`, `.samplingX`/`Y`, `.reciprocalSamplingX`/`Y`, `.matches(_:)` | **partial**: gpts and extent only (no `sampling=` input; two dimensions; endpoint False); gpts restricted to f · 2^m with f in {1, 3, 5, 15} and m ≥ 3 (abTEM takes any); the reciprocal sampling is abTEM's expression (bit-equal, slice-2 refuter B1) |
| `core.grid.spatial_frequencies`, NumPy `fftfreq` | `Grid.fftfreq`, `.frequenciesX`/`Y` | **done** (exact today, `SimulatoCoreTests`; D12 allows a tolerance) |
| `BaseWaves.angular_sampling`, `.cutoff_angles`, `.antialias_cutoff_gpts`, `.antialias_valid_gpts`, `.rectangle_cutoff_angles`, `.full_cutoff_angles`, `.cutoff_frequencies` | `DiffractionIntensity.angularSampling`, `.cutoffAngles`, `Antialias.cutoffPatternSize(on:)` | **partial**: angular sampling, cutoff angles and the odd cutoff size; the valid size only inside `PixelatedDetector(.valid)`; the rectangle and full cutoff angles and the cutoff frequencies not offered |
| `abtem.config` (dask, device, precision, antialias settings), `abtem.axes` (axis metadata) | | **absent**: not planned: global configuration and ensemble-axis metadata belong to abTEM's lazy array objects; the kit takes explicit arguments |
| `abtem.distributions` (parameter distributions for ensembles) | | **absent**: not planned: a caller loops over a parameter and averages (focal series, tilt series) |

## Potentials and slicing

| abTEM | SimulatoKit | State |
|---|---|---|
| `Potential(atoms, gpts, slice_thickness, parametrization="lobato", projection="infinite")`, the lazy potential, and its `device` | `Potential(structure:nx:ny:sliceThickness:parameterization:)`, `.build() -> ProjectedPotential` (registration-performance.md section 1) | **done**: the grid, the slice thicknesses, the wrapped atoms with their slice and the parameterisation, without the arrays; `build()` equals the former eager construction (asserted exactly today; D12 allows a tolerance); on `.metal` the potential is built on the GPU (V within 1e−5 × max V of the CPU's, section 3 there) |
| `ase.Atoms` as the atoms argument (cell, positions, numbers, `pbc`) | `AtomicStructure(box:atoms:)`, `.repeated(_:_:_:)` | **partial**: orthogonal (diagonal) boxes, atomic number and Cartesian position only; periodic on all axes (abTEM honours ASE's `pbc`: amendment 3, item 5); repeats in ASE's order; no ASE or CIF reader in the library; a local type, marked for replacement by DiffractoKit's cell types |
| `Potential(slice_thickness=float or sequence)`, `SliceIndexedAtoms` | `sliceThickness:`, `sliceThicknesses:`, `Slicing.thicknesses`, `.validated`, `.wrapped`, `.sliceIndices` | **done**: ⌈Lz / s⌉ equal slices, or a list that must sum to Lz (NumPy `isclose`), periodic wrap, top-face rule; no default of 1 Å (the thickness is required); an atom in the gap of a short slice list is kept, abTEM drops it (amendment 3, item 3) |
| `Potential(sampling=)` | | **absent**: planned: SimulatoMultislice (the grid from a sampling) |
| `Potential(exit_planes=)` | `ExitPlanes` on `Multislice.run` | **done**: `.every(m)` includes the entrance plane as abTEM does; `.afterSlices` always adds the last slice, abTEM's tuple does not (amendment 3, item 6) |
| `Potential(projection="finite")` (`QuadratureProjectionIntegrals`), `SlicedAtoms(xy_padding, z_padding)` | | **absent**: planned: SimulatoMultislice (finite projections for thick slices) |
| `Potential(periodic=False)` (cut cell with margins) | | **absent**: planned: SimulatoMultislice |
| `Potential(plane, origin, box)` (rotated and cut cells, `orthogonalize_cell`) | | **absent**: planned: after DiffractoKit's cell types (`plane="xy"`, `origin` 0 and the cell's own box only) |
| `Potential(integrator=)`, `ScatteringFactorProjectionIntegrals`, `GaussianProjectionIntegrals`, `ProjectionIntegralTable` | | **absent**: not planned: the infinite scattering-factor projection is the only integrator |
| `Potential.build(first_slice, last_slice, max_batch, lazy)`, `.generate_slices`, iteration | `Potential.build()` | **partial**: the whole potential only; `first_slice`, `last_slice`, `max_batch`, `lazy` (dask) and slice generation not offered |
| `PotentialArray(array, slice_thickness, extent, sampling, exit_planes)` (the eager potential) | `ProjectedPotential` (`grid`, `slices`, `sliceOfAtom`, `init(grid:slices:)`) | **partial**: built by `Potential.build()` (done); an array from outside cannot be supplied: `PotentialSlice` has no public initialiser; no ensemble axes or stored exit planes |
| `PotentialArray.num_slices`, `.slice_thickness`, `.thickness`, `.slice_limits`, `.extent`, `.sampling` | `ProjectedPotential.slices.count`, `PotentialSlice.thickness`, `.thickness`, `.grid` | **partial**: `slice_limits` not offered |
| `PotentialArray.transmission_function(energy)` | `Multislice.transmissionFunctions(_:beam:)` | **done**: the bandlimited t = IFFT(FFT(exp(iσV)) · A) |
| `TransmissionFunction.transmit(waves, conjugate)` | | **absent**: planned: SimulatoMultislice (the transmission step on its own; it runs inside `Multislice.run`) |
| `PotentialArray.project()`, `.to_images()` | | **absent**: planned: SimulatoMultislice (the summed slices as a `RealField`) |
| `PotentialArray.tile`, `CrystalPotential` (tiled potentials) | | **absent**: planned: SimulatoMultislice |
| `PotentialArray.depth_profile`, `.show`, `.show_depth_profile` | | **absent**: not planned: plotting (the kits have no UI) |
| `crystal_slice_thicknesses(atoms, tolerance)` | | **absent**: planned: SimulatoMultislice (`Slicing`) |
| `ChargeDensityPotential`, `GPAWPotential` | | **absent**: not planned: needs DFT inputs; outside the field's common use of the suite |
| `orthogonalize_cell`, `standardize_cell`, the helpers of `atoms.py` (`rotate_atoms`, `cut_cell`, `pad_atoms`, `shrink_cell`, ...) | | **absent**: planned: with DiffractoKit's cells |

## Parametrizations

| abTEM | SimulatoKit | State |
|---|---|---|
| `LobatoParametrization` (`projected_scattering_factor`, `scaled_parameters`) | `LobatoParameterization` (`projectedScatteringFactor(atomicNumber:k2:)`, static `scaledParameters(atomicNumber:)`, `parameters(atomicNumber:)`, `supports(atomicNumber:)`) | **partial**: neutral atoms (`charge` not offered), Z = 1 to 103, the shipped table only (`parameters=` custom JSON and `sigmas=` not offered); `.abTEMFloat32` bit-close, `.exact` (the default) registered at 1e−7 |
| `LobatoParametrization.potential`, `.projected_potential`, `.scattering_factor`, `.charge`, `.x_ray_scattering_factor`, `.finite_projected_potential`, `.line_profiles` | | **absent**: planned: SimulatoMultislice (the closed forms; `projected_potential` is the far-field truth's reference in the tests, evaluated in Python) |
| `KirklandParametrization`, `PengParametrization`, `EwaldParametrization` | | **absent**: planned: SimulatoMultislice |
| `Parametrization.sigmas`, `sigma_to_B`, `B_to_sigma`, `validate_sigmas` (Debye–Waller tables) | | **absent**: planned: SimulatoMultislice (σ is given by the caller per element today) |
| `Parametrization(parameters)`, `.to_json`, `.from_json`, `.fit` | | **absent**: not planned: custom tables and fitting are research tools |

## Waves and multislice

| abTEM | SimulatoKit | State |
|---|---|---|
| `PlaneWave(energy, extent, gpts, sampling, normalize=False, device)`, `.build()` | `Waves.planeWave(beam:grid:)` | **partial**: `normalize=True` not offered; the grid is given as gpts and extent; the backend is chosen at `Multislice.run` |
| `PlaneWave(tilt=)`, `Probe(tilt=)`, `Waves.base_tilt` | | **absent**: planned: SimulatoMultislice |
| `Waves(array, energy, extent, sampling, reciprocal_space, ensemble_axes_metadata, metadata)` | `Waves(field:grid:beam:normalization:)` | **partial**: one two-dimensional wave in real space; no ensemble axes, reciprocal-space storage or metadata; `normalization` is `.values` or `.reciprocalSpace` |
| `Probe(semiangle_cutoff, energy, extent, gpts, sampling, soft=True, aperture, aberrations)` | `Probe(beam:grid:aperture:aberrations:)` | **partial**: no `tilt`, `scan_positions`, `metadata`, distributions of parameters, or aberration coefficients as keywords (set on `Aberrations`); a probe is rebuilt to change a parameter |
| `Probe.build(scan=)` | `Probe.build(at:)` | **partial**: one position per call (a scan: `STEM.scan`); `scan=` as a grid of positions not offered |
| `Probe.scan(potential, scan, detectors)`, `Probe.multislice(potential, scan, detectors)` | `STEM.scan(probe:potentials:scan:detectors:backend:)` | **done**: slice 2; positions bit-equal, detectors averaged over the potentials; `max_batch` is chosen inside (chunks of at most 256 MB of waves), `lazy` not offered |
| `Waves.multislice(potential, detectors)`, `PlaneWave.multislice(...)` | `Multislice.run(_:through:exitPlanes:backend:)` | **partial**: waves in, waves at the exit planes out; no `detectors=` (apply a `Detector` to `DiffractionIntensity(exit)`, as `STEM.scan` does per position) |
| `multislice_and_detect(waves, potential, algorithm=FourierMultislice(order=1))` | `Multislice.run(_:through:exitPlanes:backend:)` | **done**: `.cpu` float64 and `.metal` float32 (Backends) |
| `FourierMultislice(order=2)`, `conjugate`, `transpose` | | **absent**: planned: SimulatoMultislice |
| `RealSpaceMultislice` (finite differences, backscattering) | | **absent**: not planned for now (a research option) |
| `FresnelPropagator.get_array(waves, thickness, order=1)` | `Multislice.propagator(grid:beam:thickness:)` | **done**: order 1, bandlimited; pinned through the exit waves |
| `FresnelPropagator.propagate(waves, thickness, in_place, order)` | | **absent**: planned: SimulatoMultislice (free propagation of a wave over Δz; order 2 with `FourierMultislice(order=2)`) |
| `AntialiasAperture`, `antialias_aperture(gpts, sampling)`, `bandlimit` | `Antialias.aperture(on:)`, `.cutoff`, `.taper` | **done**: abTEM's defaults 0.6666666 and 0.01 (not configurable) |
| `Waves.intensity()` | `Waves.intensity()` | **done** |
| `Waves.diffraction_patterns(max_angle="cutoff" or "full", fftshift=True)` | `Waves.diffractionPattern(maxAngle:)` | **partial**: `"valid"`, a numeric angle, `block_direct`, `return_complex`, `renormalize=False`, `fftshift=False` and `parity` other than "odd" not offered |
| `Waves.phase()`, `.real()`, `.imag()`, `.to_images()` | `ComplexField.real`, `.imaginary`, `.squaredMagnitude()` | **partial**: raw arrays; no phase helper and no image-returning methods |
| `Waves.apply_ctf(ctf)` (HRTEM image formation) | `Waves.applying(_:)` | **done**: slice 3; the `**kwargs` shorthand that builds a CTF is not offered |
| `Waves.normalize`, `.tile`, `.downsample`, `.convolve`, `.phase_shift`, `.ensure_reciprocal_space`, `.ensure_real_space` | | **absent**: planned as needed: SimulatoCore |
| `Waves.scan`, `.depth_profile`, `.show_depth_profile`, `.show`, `Probe.profiles`, `.ctf`, `.show` | | **absent**: not planned: plotting and lazy-array helpers |
| `Waves.transition_potential_multislice`, `Probe.transition_potential_scan` | | **absent**: not planned: inelastic scattering (see Other abTEM capabilities) |

## Probe and transfer

| abTEM | SimulatoKit | State |
|---|---|---|
| `Aperture(semiangle_cutoff, soft=True)`, `soft_aperture`, `hard_aperture` | `Aperture(semiangleCutoff:soft:)`, `.values(on:wavelength:)` | **done** |
| `Aberrations` coefficients C10 to C56 and φ (25 symbols), χ(α, φ) | `Aberrations` (`C10` ... `phi56`), `.chi(alpha:phi:wavelength:)` | **done**: slice 1 section 3.8, orders 1 to 5 |
| `Aberrations` aliases: `defocus` = −C10, `Cs` = C30, and 23 more (`C5`, `astigmatism`, `astigmatism_angle`, `coma`, `trefoil`, `quadrafoil`, `pentafoil`, `hexafoil`, `astigmatism3`, `astigmatism5`, `coma4`, `trefoil4`, `quadrafoil5` and their `_angle`) | `Aberrations.defocus`, `.Cs` | **partial**: `defocus` and `Cs` only; the other 23 by symbol |
| `Aberrations(aberration_coefficients=dict, **kwargs)`, `set_aberrations`, `.aberration_coefficients` | | **absent**: not planned: the coefficients are stored properties, set directly |
| `CTF(semiangle_cutoff, soft, focal_spread, angular_spread, aberration_coefficients, energy, extent, gpts, sampling, flip_phase, wiener_snr)` | `ContrastTransferFunction(aberrations:semiangleCutoff:softAperture:focalSpread:angularSpread:)`, `.kernel(on:wavelength:)` | **partial**: slice 3; `flip_phase` and `wiener_snr` not offered; energy and grid are given at `kernel(on:wavelength:)`, not stored; negative spreads refused (abTEM overflows) |
| `SpatialEnvelope(angular_spread, aberration_coefficients)`, `TemporalEnvelope(focal_spread)` (partial coherence) | `ContrastTransferFunction.spatialEnvelope`, `.temporalEnvelope` | **done**: slice 3, values at (α, φ) and inside the kernel; the quasi-coherent approximation, as abTEM |
| `transfer.scherzer_defocus`, `point_resolution`, `nyquist_sampling` (functions) | `ContrastTransferFunction.scherzerDefocus`, `.pointResolution`, `.nyquistSampling` | **done**: slice 3 |
| `CTF.scherzer_defocus`, `.point_resolution`, `.crossover_angle` (from the CTF's own Cs and energy), `BaseAperture.nyquist_sampling` | | **absent**: planned: SimulatoProbe (instance forms; the statics take Cs and energy) |
| `CTF.to_point_spread_functions`, `.profiles`, `.to_diffraction_patterns` | | **absent**: planned: SimulatoProbe |
| `CTF.show`, `BaseTransferFunction.show` | | **absent**: not planned: plotting |
| `Bullseye`, `Vortex`, `AnnularAperture`, `Zernike`, `RadialPhasePlate` | | **absent**: planned as needed: SimulatoProbe |
| `transfer.polar2cartesian`, `cartesian2polar` | | **absent**: planned: SimulatoProbe |

## Frozen phonons

| abTEM | SimulatoKit | State |
|---|---|---|
| `FrozenPhonons(atoms, num_configs, sigmas, directions, ensemble_mean, seed)` | `FrozenPhonons(structure:sigmas:configurations:seed:directions:)`, `FrozenPhonons(structure:displacements:)` | **partial**: see the five rows below |
| `sigmas` | `sigmas: [Int: Double]` (Å, per atomic number) | **partial**: a σ for every element of the structure, as abTEM requires; a scalar for all atoms, symbols as keys, per-atom arrays and anisotropic σ (a 3-vector) not offered |
| `num_configs` | `configurations:`, `FrozenPhonons.count` | **done** |
| `seed` | `seed: UInt64`, `FrozenPhonons.configurationSeed(_:_:)` | **partial**: abTEM's displacements reproduced only when given (`displacements:`, bit-equal); NumPy's PCG64 normals are not ported, the kit's generator (`GaussianGenerator`: SplitMix64, xoshiro256**, Box–Muller) agrees in distribution (KS tests); a seed is required (no None) |
| `directions="xyz"` | `FrozenPhonons.Directions` (`.x`, `.y`, `.z`, `.all`) | **done** |
| `ensemble_mean=True` | averaging inside `STEM.scan` | **partial**: the mean only; `ensemble_mean=False` (one measurement per configuration) not offered: planned: SimulatoDetector (a caller can scan each potential alone) |
| `FrozenPhonons.atoms`, `.randomize`, `.num_configs`, `Potential(frozen_phonons)` | `FrozenPhonons.configuration(_:)`, `.displacements`, `.potentials(nx:ny:sliceThickness:parameterization:)` (lazy: `lazyPotentials`, registration-performance.md section 1) | **done**: slice 2; one potential per configuration, averaged by `STEM.scan` |
| `FrozenPhonons.to_atoms_ensemble`, `AtomsEnsemble(trajectory)`, `DummyFrozenPhonons` | | **absent**: planned: SimulatoMultislice (a trajectory of structures, for example from molecular dynamics) |

## Scans

| abTEM | SimulatoKit | State |
|---|---|---|
| `GridScan(start, end, gpts, endpoint=False)` | `GridScan(start:end:nx:ny:endpoint:)`, `.positions`, `GridScan.linspace` | **done**: slice 2; positions bit-equal to abTEM's (NumPy's `linspace`, x the first index); a scan of extent ≤ 0 on both axes refused, as abTEM |
| `GridScan(endpoint=True)` | `endpoint: true` | **partial**: implemented as `linspace(endpoint: true)`; the pins use `endpoint=False`, so it is not pinned against abTEM |
| `GridScan(sampling=)` | | **absent**: planned: SimulatoProbe (gpts from a sampling) |
| `GridScan(potential=, fractional=, end=None)` (extent from the potential, fractional coordinates, `Atom` start and end) | | **absent**: planned: SimulatoProbe (`start` and `end` are required, in Å) |
| `LineScan(start, end, gpts, sampling, endpoint=True, fractional, potential)`, `.at_position`, `.add_margin` | | **absent**: planned: SimulatoProbe |
| `CustomScan(positions)` | | **absent**: planned: SimulatoProbe (`Probe.build(at:)` takes any position today; `STEM.scan` takes a `GridScan`) |
| `BaseScan.get_positions`, `.limits`, `.num_positions`, `.shape` | `GridScan.positions`, `.nx`, `.ny` | **partial**: positions and the shape; `limits` not offered |
| `GridScan.add_to_plot`, `LineScan.add_to_plot` | | **absent**: not planned: plotting |

## Detectors

| abTEM | SimulatoKit | State |
|---|---|---|
| `AnnularDetector(inner, outer, offset)` | `AnnularDetector(inner:outer:)` | **partial**: slice 2; `offset` not offered; a negative inner is refused (abTEM squares it); `to_cpu`, `url` do not apply |
| `PixelatedDetector(max_angle, resample, reciprocal_space)` | `PixelatedDetector(maxAngle:)` | **partial**: "cutoff", "valid", "full"; numeric angles, `resample` and real-space detection not offered |
| `SegmentedDetector(nbins_radial, nbins_azimuthal, inner, outer, rotation, offset)` | `SegmentedDetector(radialBins:azimuthalBins:inner:outer:rotation:)` | **partial**: summed in float64 where abTEM sums in float32 (slice-2 amendment 1); inner = outer refused; `offset` not offered |
| `FlexibleAnnularDetector(step_size, inner, outer)` | | **absent**: planned: SimulatoDetector |
| `WavesDetector(gpts)` | | **absent**: planned: SimulatoDetector (the complex exit wave per position; `Multislice.run` returns it for one probe today) |
| `BaseDetector.detect(waves)`, `.angular_limits` | `Detector.detect(_:)` on a `DiffractionIntensity`, `.shape(for:beam:)` | **done**: slice 2 |
| `AnnularDetector.get_detector_region`, `get_detector_regions`, `.show` | | **absent**: not planned: plotting |

## Measurements

| abTEM | SimulatoKit | State |
|---|---|---|
| `Images(array, sampling)` | `RealField` (`nx`, `ny`, `values`, `samplingX`, `samplingY`, `sum`, `max`), `ScanMeasurement.image(samplingX:samplingY:)` | **partial**: a two-dimensional array with its sampling; no ensemble axes, metadata or units |
| `DiffractionPatterns(array, sampling, fftshift)` | `RealField` from `Waves.diffractionPattern`, `DiffractionIntensity` (unshifted, with `angularSampling`) | **partial**: the pattern with its sampling in 1/Å; no angular sampling or units stored on the `RealField` |
| `PolarMeasurements` (`integrate`, `integrate_radial`, `differentials`, `to_diffraction_patterns`) | `ScanMeasurement` of a `SegmentedDetector` (values over radial × azimuthal bins) | **absent**: planned: SimulatoDetector (the operations; the bins exist as raw values) |
| `RealSpaceLineProfiles`, `ReciprocalSpaceLineProfiles`, `_BaseMeasurement1D.width`, `.interpolate` | | **absent**: planned: SimulatoDetector (with `LineScan`) |
| `MeasurementsEnsemble`, `BaseMeasurements.to_measurement_ensemble`, `normalize_ensemble` | | **absent**: not planned: Swift arrays; a caller keeps the configurations apart |
| `IndexedDiffractionPatterns`, `DiffractionPatterns.index_diffraction_spots` | | **absent**: planned: with DiffractoKit's indexing |
| `BaseMeasurements.intensity()` | `Waves.intensity()`, `ComplexField.squaredMagnitude()` | **done** |
| `BaseMeasurements.real()`, `.imag()` | `ComplexField.real`, `.imaginary` | **partial**: raw arrays, not measurements |
| `BaseMeasurements.phase()`, `.abs()` | | **absent**: planned: SimulatoCore |
| `BaseMeasurements.reduce_ensemble()` | the mean in `STEM.scan` | **partial**: the mean over configurations only |
| `BaseMeasurements.relative_difference(other, min_relative_tol)` | | **absent**: planned: SimulatoCore |
| `poisson_noise(dose_per_area, total_dose, samples, seed)` (`BaseMeasurements`, `DiffractionPatterns`, `PolarMeasurements`) | | **absent**: planned: SimulatoDetector (with detector noise; a simulated image is not yet a micrograph) |
| `DiffractionPatterns.integrate_radial(inner, outer, offset)` | `AnnularDetector.detect(_:)` | **partial**: only as a detector on an exit wave's `DiffractionIntensity`; not a method of a stored pattern; `offset` not offered |
| `DiffractionPatterns.polar_binning(nbins_radial, nbins_azimuthal, inner, outer, rotation, offset)` | `SegmentedDetector.detect(_:)` | **partial**: as the row above |
| `DiffractionPatterns.radial_binning(step_size, inner, outer)`, `.azimuthal_average` | | **absent**: planned: SimulatoDetector |
| `DiffractionPatterns.center_of_mass(units)`, `.integrated_center_of_mass()` | | **absent**: planned: SimulatoDetector |
| `DiffractionPatterns.crop(max_angle, max_frequency, gpts)` | `Waves.diffractionPattern(maxAngle: .cutoff / .full)`, `PixelatedDetector(maxAngle:)` | **partial**: the crop is chosen when the pattern is made; not a method of a stored pattern |
| `DiffractionPatterns.interpolate(sampling, gpts)`, `Images.interpolate(sampling, gpts, method, boundary, order, normalization, cval)` | | **absent**: planned: SimulatoCore |
| `DiffractionPatterns.bandlimit`, `.block_direct`, `.gaussian_source_size`, `.tile_scan` | | **absent**: planned: SimulatoDetector |
| `Images.crop(extent, offset, centered)` | | **absent**: planned: SimulatoCore |
| `_BaseMeasurement2D.gaussian_filter(sigma, boundary, cval)` | | **absent**: planned: SimulatoCore |
| `interpolate_line`, `interpolate_line_at_position`, `integrate_disc`, `Images.integrate_gradient`, `Images.tile` | | **absent**: planned: SimulatoCore |
| `Images.scan_noise(dwell_time, flyback_time, rms_power, ...)`, `Images.diffractograms()` | | **absent**: not planned for now: scan-noise modelling; diffractograms are the field of an image-analysis kit |
| `BaseMeasurements.show`, `to_zarr`, `to_tiff`, `to_hyperspy`, `to_data_array`, `to_quantem` | | **absent**: not planned: plotting and file formats |

## Other abTEM capabilities

| abTEM | SimulatoKit | State |
|---|---|---|
| `SMatrix`, `SMatrixArray` (PRISM) | | **absent**: not planned in 2026 (a later acceleration lane) |
| `BlochWaves`, `StructureFactor` | | **absent**: not planned for now (DiffractoKit owns structure factors; Bloch waves a possible later slice) |
| `inelastic.core_loss`, `inelastic.plasmons`, `transition_potential_*` | | **absent**: not planned: SpectroKit's domain if ever |
| `mtf.py` (detector MTF), `noise.py` | | **absent**: planned: SimulatoDetector (with the noise rows above) |
| `reconstruct.py` (ptychography), `mcf.py`, `magnetism`, `finite_difference.py` | | **absent**: not planned: research modules outside the kit's scope |
| `array.stack`, `concatenate`, `from_zarr`, dask laziness | | **absent**: not planned: Swift arrays; laziness is FormatoKit's |
| `visualize`, `show_atoms` | | **absent**: not planned: the kits have no UI |

## Backends

abTEM runs on NumPy, or on CuPy with `device="gpu"`. SimulatoKit has two backends behind one API (`Backend`, passed to
`Multislice.run`, `Multislice.runBatch` and `STEM.scan`; a `Potential` on `.metal` builds on the GPU, registration-performance.md section 1):

| | Precision and engine | Role | Agreement |
|---|---|---|---|
| `.cpu` | float64, Accelerate; a batch of positions on at most 8 threads, equal to sequential runs (exact today; D12 allows a tolerance) | the reference | against abTEM 1.0.10 (float64 configuration) at the pins' bounds, `docs/registration-slice-1.md` section 4 (V 1e−11 × max V, exit waves, intensities and patterns 1e−10 × max with `.abTEMFloat32`; 1e−6 with `.exact`), slice 2 (annular and pixelated 1e−10, segmented 1e−6) and slice 3 (kernels 1e−12, images 1e−10) |
| `.metal` | float32, MPSGraph FFT | abTEM's `device="gpu"` counterpart | against `.cpu`, registered bounds: `docs/registration-slice-1.md` section 4 (exit wave 1e−4, intensity and pattern 2e−4 of the maximum; Σ\|ψ\|² 1.5e−6 per slice, amendment 2), slice 2 (scans 2e−4), and `docs/registration-performance.md` section 3 (V built on the GPU 1e−5 × max V, GPU probe 1e−5 × max \|ψ\|, grids to 1024²) |

The Metal bounds hold to 100 slices (slice-1 amendment 3, item 7); MPSGraph's float32 FFT is not exactly unitary, a known
limit (amendment 2; `docs/status.md`). The CTF (`Waves.applying`) runs on the CPU only. A Metal result beyond a bound is a
failure of the Metal path, never a tolerance to widen. `MetalAccuracy.registered` carries these bounds as a value.

Execution choices are parameters (`ExecutionOptions`, passed as `options:`; 0.2.0, `docs/registration-generality.md`):
threads, Metal batch memory and shapes, the CPU scan chunk, the renormalisation, the executable cache. abTEM's nearest
counterparts are its `max_batch` and dask settings, which are not counted as missing options above.

## Not abTEM's, SimulatoKit's own

- The Metal backend (`Backend.metal`) and the batched scan (`Multislice.runBatch`), registered against the CPU path; MPSGraph's float32 FFT gain is a known limit
  (registration amendment 2).
- The CPU batch on up to `ExecutionOptions.maximumCPUThreads` (8 by default) threads, equal to sequential runs (exact today; D12 allows a tolerance).
- The constants of the CTF closed forms as parameters (`factor:` of `scherzerDefocus`, 1.5 by default as abTEM; `divisor:` of
  `pointResolution`, 6; `oversampling:` of `nyquistSampling`, 4).
- `LobatoParameterization.PiConstant`: the exact π by default, abTEM's float32 π on request.
