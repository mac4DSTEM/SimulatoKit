# tools/metal-fft-gain: is MPSGraph's float32 FFT round trip unitary?

A diagnostic, not a test. It explains the Metal path's total-intensity bound per slice
(`MetalAccuracy.registered.totalIntensityPerSlice`, 1.5e−6) and tells you when that bound can be tightened.

What it found (2026-10-15; macOS 27, Apple M5 Pro, Xcode 27.0): MPSGraph's complex float32 FFT round trip (forward `.none`,
inverse `.size`, or both `.unitary`) is not exactly unitary. Each round trip multiplies the non-zero frequencies by a gain
of about 1 − 5e−7 (Σ|ψ|² falls by 1.4e−7 to 1.4e−6 per round trip depending on the grid and the field), the zero frequency
passes exactly, and the remainder after a fitted scalar is 1e−13 of the energy. It is deterministic and linear in the
number of round trips. One multislice step is one round trip, so the Metal path's total intensity drifts below the CPU
path's by about 8e−7 per slice at 512² without the per-step renormalisation (`ExecutionOptions.metalRenormalization`).

```
xcrun swiftc -O tools/metal-fft-gain/round-trip-drift.swift -o /tmp/drift && /tmp/drift   # Σ|ψ|² after 1, 10, 50, 100 round trips
xcrun swiftc -O tools/metal-fft-gain/round-trip-gain.swift -o /tmp/gain && /tmp/gain      # the fitted gain per input
```

Re-run both after a macOS or Xcode update. If the gain is gone (Σ|ψ|² holds to float32 rounding), the total-intensity
bound can return to 1e−5 × the maximum like the other Metal bounds. Output is printed only; logs are not kept in the repo.
