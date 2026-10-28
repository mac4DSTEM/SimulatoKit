# Contributing to SimulatoKit

Thank you for your interest. Issues and pull requests are welcome.

## Build and test

```sh
swift test -c release
```

Build in release mode: the FFT and potential loops are many times slower unoptimised. Do not use `swift build --build-tests -c release`; it fails on the current toolchain. The GPU (Metal) tests
skip themselves with a stated reason when no Metal device is present.

## How numbers are checked

- The CPU path is float64 and is compared with [abTEM](https://github.com/abtem/abtem) 1.0.10, using fixtures generated
  by the scripts in `tools/slice-*/pins.py`. A test reads the fixture and states its tolerance.
- A test's tolerance is the one the Python package's own tests use for that function; a result outside it is a bug. A change
  that moves a fixture says why in the pull request and in the CHANGELOG. If you believe abTEM is wrong, open an issue
  with the evidence.
- A change that moves a number needs a test that fails when the code is broken, and a comparison against abTEM (or a
  closed form); say which in the pull request.
- The Metal (float32) path is compared with the CPU path within the bounds in the README.

## Pull requests

- Keep a change focused, with tests that fail before and pass after.
- Update the README, `docs/parity.md` and `CHANGELOG.md` when behaviour or the public API changes.
- The library is GPL-3.0-only; by contributing you agree your work is released under it.

## Reporting issues

Please include the SimulatoKit version, macOS version and hardware, the backend (`.cpu` or `.metal`), and a minimal
snippet that shows the problem.
