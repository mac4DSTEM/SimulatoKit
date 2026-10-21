# Status

SimulatoKit 0.3.1 (2026-10-21), against abTEM 1.0.10. Live: crystal cells and the Lobato projected potential, any grid length of at least 2 (vDSP where it has
the length, Bluestein elsewhere, slower than NumPy there), multislice on `.cpu` (float64) and `.metal` (float32, held to the CPU), beam tilt, the probe with
polar aberrations C10 to C56, `GridScan`, frozen phonons, annular, pixelated, segmented and centre-of-mass detectors, `STEM.scan`, the contrast transfer function,
source size and Poisson noise. What is checked against what: `docs/parity.md` (28 present, 36 partial, 59 missing of 123 rows); timings against Python:
`docs/benchmarks.md`; what is open: `docs/open-items.md`.
