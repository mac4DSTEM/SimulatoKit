# Status

SimulatoKit 0.4.0 (2026-10-22), against abTEM 1.0.10. Live: crystal cells and the Lobato projected potential, any grid length of at least 2 (vDSP where it has
the length, Bluestein elsewhere, slower than NumPy there), multislice on `.cpu` (float64) and `.metal` (float32, held to the CPU), beam tilt, the probe with
polar aberrations C10 to C56, `GridScan`, frozen phonons, annular, pixelated, segmented and centre-of-mass detectors, `STEM.scan`, the contrast transfer function,
source size and Poisson noise, `LineScan` and line profiles, integrated centre of mass (`integrate_gradient`, iCOM), and `WavesDetector` (the exit wave at every scan position). What is checked against what: `docs/parity.md` (37 present, 36 partial, 55 missing of 128 rows); timings against Python:
`docs/benchmarks.md`; what is open: `docs/open-items.md`.

Tests are split in two. `SimulatoKitParityTests` holds the pinned abTEM comparisons, the truth tests, Metal against CPU and the benchmarks (111 tests; fixtures live in its `Fixtures/`); the other targets hold the fast unit tests (35 tests); `SimulatoKitTestSupport` holds the helpers both use. `swift test -c release --parallel --skip SimulatoKitParityTests` runs the fast half; the plain `swift test -c release --parallel` runs all 146.
