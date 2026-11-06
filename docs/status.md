# Status

SimulatoKit 0.5.0 (2026-11-06), against abTEM 1.0.10. Live: crystal cells and the Lobato projected potential (`project`, tiled; slices that follow the crystal planes), any grid length of at least 2 (vDSP where it has
the length, a mixed-radix FFT for 2, 3, 5 and 7 and a generic odd-radix pass for the other lengths, Bluestein where a large prime factor makes it the cheaper), multislice on `.cpu` (float64) and `.metal` (float32, held to the CPU), beam tilt, the probe with
polar aberrations C10 to C56, `GridScan`, frozen phonons (exit waves too), annular, pixelated, segmented (DPC: integrate, integrateRadial, differentials on a polar scan) and centre-of-mass detectors, `STEM.scan`, the contrast transfer function,
source size and Poisson noise (samples too), `LineScan` (matched to a probe: end and Nyquist sampling defaults) and line profiles, integrated centre of mass (`integrate_gradient`, iCOM), the azimuthal average of a centred pattern, Fourier-space image interpolation and `integrateDisc(interpolate:)`, and `WavesDetector` (the exit wave at every scan position). What is checked against what: `docs/parity.md` (38 present, 40 partial, 51 missing of 129 rows); timings against Python:
`docs/benchmarks.md`; what is open: the GitHub issues; known limits: `docs/parity.md`.

Since 0.4.0: 0.5.0 (2026-10-31) added the items above and the mixed-radix and odd-radix FFT passes (250² 2.2x, 500² 2.3x off vDSP). Details: `CHANGELOG.md`.

Tests are split in two. `SimulatoKitParityTests` holds the pinned abTEM comparisons, the truth tests, Metal against CPU and the benchmarks (111 tests; fixtures live in its `Fixtures/`); the other targets hold the fast unit tests (35 tests); `SimulatoKitTestSupport` holds the helpers both use. `swift test -c release --parallel --skip SimulatoKitParityTests` runs the fast half; the plain `swift test -c release --parallel` runs all 146.
