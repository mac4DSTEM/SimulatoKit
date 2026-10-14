# CLAUDE.md — start here

SimulatoKit: the macEM suite's electron-scattering simulation library, a Swift port of abTEM (projected potentials,
multislice, probes; detectors and frozen phonons next), checked against abTEM 1.0.10 in float64, the truth source for the
other kits' tests. GPL-3.0, crediting abTEM. Modules: `SimulatoCore` (beam, grid, fields, waves, the vDSP float64 FFT,
the antialias aperture, the local structure type), `SimulatoMultislice` (slice 1: Lobato, slicing, potential, multislice on
`.cpu` and `.metal`), `SimulatoProbe` (slice 1), `SimulatoDetector` (placeholder, slice 2). Rules only here; state in
`docs/status.md`.

## Read first

1. The suite's `../CLAUDE.md` (board, tiers, how we work) and `../ROADMAP.md`, "Kits" (the vision, the Apple-hardware rule,
   phase 5).
2. `docs/registration-slice-1.md`: the definitions, the pins, the truth, the gates, the two amendments.
3. `docs/status.md`: what is live, the hand-off, what is unverified. `docs/parity.md`: abTEM's API, ported or not.
4. `../SpectroKit/CLAUDE.md`: the model for a kit repo (this file follows it).

## Hard rules

- **The apps and the other kits are read-only from here.** Never build, test or write in them.
- **Registration before code** for anything tier 1: `docs/registration-<slice>.md` with definitions, pins, planted truth,
  gates and planned mutants, committed before the first Swift line; a change after code starts is a dated amendment.
- **The gate is `tools/test.sh`** (`swift test -c release "$@"`). Never run `swift build --build-tests -c release`.
- **Green stamp**: a full, passing `tools/test.sh` on a clean tracked tree writes `.build/last-green`; `tools/last-green.sh`
  prints it and says whether HEAD equals it.
- **Hook**: `.claude/settings.json` runs `tools/hooks/no-sleep-poll.sh` before every Bash call (no sleep of 30 s or more, no
  while/until loop with a sleep). Wait on notifications or Monitor. Proof: `tools/hooks/test-no-sleep-poll.sh`.
- **A reference for every number.** abTEM 1.0.10 (tag `v1.0.10` = `164e644`), ASE 3.29.0, NumPy 2.5.4, SciPy 1.18.1, numba
  0.68.0, Python 3.14.7, in the gitignored `References/` (`abtem-venv`, `abtem-164e644`, `abtem-venv-freeze.txt`). The
  oracle runs as `References/abtem-venv/bin/python -I -B` with abTEM's config `precision: float64`, `fft: numpy`,
  `dask.lazy: false`. To rebuild `References/` on a new machine: `git clone --branch v1.0.10 https://github.com/abtem/abtem
  References/abtem-164e644`, `python3 -m venv References/abtem-venv`, `References/abtem-venv/bin/python -m pip install
  abtem==1.0.10`, then compare with the freeze file. An unpinned copy is never a reference. Agreement with abTEM is
  agreement, not truth: truth is a closed form, a symmetry, a conservation law or a convergence (registration section 5).
- **Pins are generated before the Swift**, by `tools/slice-1/pins.py`, which also re-derives every quantity from the
  registration's text; a difference is a failure of the port, never a tolerance to widen.
- **The π choice.** abTEM evaluates π in the Lobato factor as float32; `LobatoParameterization(piConstant: .exact)` is the
  default, `.abTEMFloat32` reproduces abTEM bit-closely. Tests against the pins use `.abTEMFloat32`.
- **Every accelerated path keeps the CPU path** as its reference with a registered tolerance. The Metal path's known limit:
  MPSGraph's float32 FFT round trip loses about 5e−7 of Σ|ψ|² per slice (amendment 2, `tools/metal-fft-gain/`).
- **`AtomicStructure` is marked for replacement** by DiffractoKit's `CrystalStructure` at DiffractoKit's first tag; the
  example's CIF reader likewise by DiffractoKit's `CIFImport`.
- **Public API only, no app types**; error messages never name an app.
- **Rigour follows risk** (`../CLAUDE.md`): tier in every commit; tier 1 is registration, every new test broken once
  (`tools/slice-1/mutants.py`, one mutant at a time on a clone), one Gate B refuter on the frozen commit, Gate D when a cause
  is unproven, a fixture.
- **Heavy work**: at most 8 jobs; the 100-slice 512² Metal test is the heaviest test (about 6 s).
- **Names**: domain + `Kit`, modules with the prefix: `SimulatoCore`, `SimulatoMultislice`, `SimulatoProbe`,
  `SimulatoDetector`.
- **Pushing, tagging and creating the remote are the owner's.** One session per repo; `main` only; commit by explicit paths
  at every green checkpoint. At each tag, `docs/consumers-x.y.z.md` says what each consumer changes.
- **The owner's data** lives in the suite's `data/` folder (`$MACEM_DATA`, default `../data`): never copied here; a manifest
  line in `../data/README.md` per file used.
- **Suite board**: this repo's row is `lanes/simulatokit`; how a session opens and closes it is in `../CLAUDE.md`.
