# SimulatoKit

Multislice TEM/STEM simulation: potentials, frozen phonons, scans, detectors, CTF; a port of abTEM. Pure Swift, no UI, standalone; the suite's rules are in `../CLAUDE.md`. Start with `README.md`,
`docs/parity.md` and the task.

## Reference
abTEM 1.0.10 in float64 (fixtures by `tools/slice-*/pins.py`), Lobato-Van Dyck (2014) parameters. CPU float64 on Accelerate; `.metal` float32 on MPSGraph FFT, checked against `.cpu` at abTEM's float32 tolerance. Fixtures are generated from it (`tools/*/pins.py`) and committed under `Tests/*/Fixtures/`; a difference inside
the reference's own tolerance is agreement, outside it a bug. Readers, I/O, indices and metadata exact.

## Build and test
```sh
swift build
swift test -c release          # the gate; Debug is many times slower on the numeric loops
tools/consumer-check/run.sh    # a separate package that compiles against the public API
```
Never `swift build --build-tests -c release` (fails on this toolchain). Real-data tests skip when the data is absent.

## Modules
`SimulatoCore` (cells, potentials, waves), `SimulatoMultislice` (`Multislice.run`, `runBatch`), `SimulatoProbe` (probes, `STEM.scan`, CTF), `SimulatoDetector` (detectors).

## Rules of the kit
- Public API only; nothing names an app or an app type; error messages never say "mac4DSTEM" or "macTEM".
- Rows and columns, not x and y. Every accelerated path keeps the CPU float64 path as its reference (D16).
- A new or changed number gets a test against the reference at its upstream tolerance, and a row in `docs/parity.md`.
- `ExecutionOptions` carries the backend and thread choices; a batch of positions runs on at most 8 threads.

## Closeout
`swift test -c release` green; `CHANGELOG.md` entry with a Migration paragraph where a consumer must change something;
version bumped in the entry; one commit on `main` whose last paragraph is the hand-off. Paul tags with
`../macEM-suite/tools/ship.sh` (tag = CHANGELOG top version).
