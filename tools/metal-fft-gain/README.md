# tools/metal-fft-gain: Gate D of 2026-10-15 (registration amendment 2)

MPSGraph's complex float32 FFT round trip (forward `.none`, inverse `.size`, or both `.unitary`) is not exactly unitary on
this machine (macOS 27, Apple M5 Pro, Xcode 27.0): each round trip multiplies the non-zero frequencies by a gain of about
1 − 5e−7 (Σ|ψ|² falls by 1.4e−7 to 1.4e−6 per round trip depending on the grid and the field), the zero frequency passes
exactly, and the remainder after a fitted scalar is 1e−13 of the energy. It is deterministic and linear in the number of
round trips. One multislice step is one round trip, so the Metal path's total intensity drifts below the CPU path's by
about 8e−7 per slice at 512².

```
xcrun swiftc -O tools/metal-fft-gain/round-trip-drift.swift -o /tmp/drift && /tmp/drift   # Σ|ψ|² after 1, 10, 50, 100 round trips
xcrun swiftc -O tools/metal-fft-gain/round-trip-gain.swift -o /tmp/gain && /tmp/gain      # the fitted gain per input
```

The logs of the 2026-10-15 runs are beside the scripts. Re-run them after an OS update: if the gain is gone, the amended
bound can return to the registered 1e−5.
