# CLAUDE.md — start here

SimulatoKit: the macEM suite's electron-scattering simulation library, a Swift port of abTEM (projected potentials,
multislice, probes, frozen phonons, scans, detectors and the CTF for HRTEM), checked against abTEM 1.0.10 in float64, the truth source for the
other kits' tests. GPL-3.0, crediting abTEM. Modules: `SimulatoCore` (beam, grid, fields, waves, the vDSP float64 FFT,
the antialias aperture, the local structure type), `SimulatoMultislice` (slice 1: Lobato, slicing, potential, multislice on
`.cpu` and `.metal`; slice 2: frozen phonons, batched propagation; the performance lane: `Potential`, the GPU engine
in `MetalEngine.swift`), `SimulatoProbe` (slice 1; slice 2: `GridScan`; slice 3: `ContrastTransferFunction`),
`SimulatoDetector` (slice 2: annular, pixelated, segmented detectors and the STEM driver). Rules only here; state in
`docs/status.md`.

## Read first

1. The suite's `../CLAUDE.md` (board, tiers, how we work) and `../ROADMAP.md`, "Kits" (the vision, the Apple-hardware rule,
   phase 5).
2. `docs/registration-slice-1.md`, `-2.md` and `-3.md`: the definitions, the pins, the truth, the gates,
   the amendments; `docs/registration-performance.md` for the Metal path on the GPU.
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
  while/until loop with a sleep; no `killall` and no `pkill -f`: kill by pid, other sessions share the machine). Wait on
  notifications or Monitor. Proof: `tools/hooks/test-no-sleep-poll.sh`.
- **A reference for every number.** abTEM 1.0.10 (tag `v1.0.10` = `164e644`), ASE 3.29.0, NumPy 2.5.4, SciPy 1.18.1, numba
  0.68.0, Python 3.14.7, in the gitignored `References/` (`abtem-venv`, `abtem-164e644`, `abtem-venv-freeze.txt`). The
  oracle runs as `References/abtem-venv/bin/python -I -B` with abTEM's config `precision: float64`, `fft: numpy`,
  `dask.lazy: false`. To rebuild `References/` on a new machine: `git clone --branch v1.0.10 https://github.com/abtem/abtem
  References/abtem-164e644`, `python3 -m venv References/abtem-venv`, `References/abtem-venv/bin/python -m pip install
  abtem==1.0.10`, then compare with the freeze file. An unpinned copy is never a reference. Agreement with abTEM is
  agreement, not truth: truth is a closed form, a symmetry, a conservation law or a convergence (registration section 5).
- **Pins are generated before the Swift**, by `tools/slice-<n>/pins.py`, which also re-derives every quantity from the
  registration's text; a difference is a failure of the port, never a tolerance to widen.
- **The π choice.** abTEM evaluates π in the Lobato factor as float32; `LobatoParameterization(piConstant: .exact)` is the
  default, `.abTEMFloat32` reproduces abTEM bit-closely. Tests against the pins use `.abTEMFloat32`.
- **Every accelerated path keeps the CPU path** as its reference with a registered tolerance. The Metal path's known limit:
  MPSGraph's float32 FFT round trip loses about 5e−7 of Σ|ψ|² per slice (amendment 2, `tools/metal-fft-gain/`), data-
  dependent (no gain map corrects it); the per-step renormalisation keeps it within 5e−7 (1.4e−7 on Al; performance amendments 1 and 3).
  Every Metal loop runs in an autorelease pool (without one a 1024² scan held 73 GB).
- **`AtomicStructure` is marked for replacement** by DiffractoKit's `CrystalStructure` at DiffractoKit's first tag; the
  example's CIF reader likewise by DiffractoKit's `CIFImport`.
- **Public API only, no app types**; error messages never name an app.
- **Rigour follows risk** (`../CLAUDE.md`): tier in every commit; tier 1 is registration, every new test broken once
  (`tools/slice-<n>/mutants.py`, one mutant at a time on a clone), one Gate B refuter on the frozen commit, Gate D when a cause
  is unproven, a fixture.
- **Heavy work**: at most 8 jobs; the 1024² × 100-slice test is the heaviest (about 25 s, mostly the CPU reference). The
  benchmark (`SIMULATO_BENCHMARK=1`, or `=metal` for the GPU side) is not part of the gate.
- **Frozen phonons**: abTEM's displacements come from NumPy's PCG64 normals, which the kit does not reproduce; tests against
  abTEM give the port abTEM's displacements (`FrozenPhonons(structure:displacements:)`), and the kit's `GaussianGenerator`
  is registered by its statistics only.
- **Pins carry exact doubles.** Arrays travel as float64 bits in the sidecars; a scalar parameter that the code amplifies
  (an aberration angle under a large χ) must travel as exact bits too, not as a long decimal in the JSON (Foundation's
  `JSONSerialization` can read those a few ulp off; slice-3 refuter).
- **Segmented detectors** sum in float64; abTEM sums them in float32 (slice-2 amendment 1), so they agree to ~2e−7.
- **Names**: domain + `Kit`, modules with the prefix: `SimulatoCore`, `SimulatoMultislice`, `SimulatoProbe`,
  `SimulatoDetector`.
- **Pushing, tagging and creating the remote are the owner's.** One session per repo; `main` only; commit by explicit paths
  at every green checkpoint. At each tag, `docs/consumers-x.y.z.md` says what each consumer changes.
- **The owner's data** lives in the suite's `data/` folder (`$MACEM_DATA`, default `../data`): never copied here; a manifest
  line in `../data/README.md` per file used.
- **Consumer check**: `tools/consumer-check/run.sh` (the whole public API from a separate package, then `readme.sh`, the
  README's code built verbatim) and `controls.sh` (negative controls); run both before a tag.
- **Every choice is a parameter** (0.2.0, `docs/registration-generality.md`): a new policy gets a parameter whose default
  is the validated value, and the defaults must reproduce every fixture bit for bit. No application or suite name in
  Sources, Tests, README, DocC, `Package.swift` or the example (the registrations and consumer notes keep the record).
- **DocC**: one catalogue per module in `Sources/<Module>/<Module>.docc`, one article each; links must resolve
  (`xcrun docc convert` per module with that module's symbol graph).
- **Suite board**: this repo's row is `lanes/simulatokit`; how a session opens and closes it is in `../CLAUDE.md`.
