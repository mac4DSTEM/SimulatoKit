# Consumers of SimulatoKit 0.0.1

The first tag: slices 1 to 3 (potentials, multislice on the CPU and Metal, probes, frozen phonons, scans, detectors, the
CTF) and the Metal performance lane (`docs/registration-performance.md`). No consumer pins SimulatoKit yet, so no consumer
changes anything. A consumer that pins 0.0.1 imports the products it needs (`SimulatoCore`, `SimulatoMultislice`,
`SimulatoProbe`, `SimulatoDetector`) at `exact: "0.0.1"`; `docs/consumers-0.1.0.md` is the note for the first API
meant to be used by other packages.

What moved since the untagged `main` some readers may have built against: `.metal` now builds the transmissions on the
GPU (and, for a `Potential`, the projected potential too) and renormalises each step, so `.metal` numbers move within
their registered bounds (closer to the CPU path: Σ|ψ|² per slice from about 9e−7 to at most 1.4e−7); `.cpu` numbers do
not move. Additive API: `Potential`, `PotentialSlice.init(thickness:values:)`, `FrozenPhonons.lazyPotentials`, and
`Multislice.run`, `Multislice.runBatch` and `STEM.scan` overloads that take `Potential`.
