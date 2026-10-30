# Benchmarks against abTEM

Per-function timings of SimulatoKit against abTEM 1.0.10 on the same inputs (Al [001], 200 keV, a plane wave or a
20 mrad probe, 10 slices of 2.025 Å, about 0.1 Å per pixel, infinite projection). One warm-up, then five timed runs;
each cell is the median in milliseconds with (min-max). Written by `tools/benchmarks/run.sh` and `merge.py`.

- Machine: Apple M5 Pro, macOS 27.0.1, Apple Swift version 6.4 (swiftlang-6.4.0.34.1 clang-2100.3.34.1)
- Python 3.12.14 (NumPy 2.5.4, SciPy 1.18.1), abTEM 1.0.10 (CPU, NumPy's FFT; no GPU)
- Kit CPU is float64 on Accelerate (one thread per call); Metal is float32 on MPSGraph; abTEM is shown in float64 and float32
- Run while other jobs used the machine (load average 6 to 11): the spread columns show it; treat ratios within 20% as equal
- Sizes: 64² to 256² small and typical; **1024² and above: not measured yet**

Re-run from a fresh clone of the kit, on a quiet Mac (`MEMORY_CAP` may name a memory-cap command prefix instead of `none`):

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r tools/requirements.txt
MEMORY_CAP=none tools/benchmarks/run.sh 64,128,192,256,512
.venv/bin/python tools/benchmarks/merge.py   # writes this file
```

## fft2+ifft2

| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |
|---|---|---|---|---|---|---|
| 58² | 0.145 (0.134-0.149) |  | 0.104 (0.104-0.106) | 0.081 (0.078-0.155) | 0.72x | 0.56x |
| 64² | 0.044 (0.041-0.062) |  | 0.056 (0.055-0.063) | 0.047 (0.045-0.08) | 1.27x | 1.07x |
| 122² | 0.427 (0.427-0.513) |  | 0.627 (0.601-0.687) | 0.516 (0.467-0.537) | 1.47x | 1.21x |
| 128² | 0.301 (0.299-0.301) |  | 0.266 (0.22-0.401) | 0.234 (0.215-0.329) | 0.88x | 0.78x |
| 186² | 2.3 (1.74-3.64) |  | 0.898 (0.615-1.52) | 0.804 (0.765-1.11) | 0.39x | 0.35x |
| 192² | 1.25 (1.12-1.34) |  | 0.353 (0.351-0.74) | 0.607 (0.458-0.936) | 0.28x | 0.49x |
| 250² | 3.75 (3.28-4.94) |  | 0.572 (0.511-0.684) | 0.939 (0.678-1.26) | 0.15x | 0.25x |
| 256² | 2.43 (1.92-2.95) |  | 1.04 (0.911-1.13) | 1.48 (0.999-1.71) | 0.43x | 0.61x |
| 506² | 20.6 (13.8-36.2) |  | 9.34 (7.85-10.1) | 4.49 (4.27-4.62) | 0.45x | 0.22x |
| 512² | 9.03 (8.82-13.5) |  | 7.29 (6.73-8.27) | 4.54 (4.45-5.47) | 0.81x | 0.50x |

(ratios are abTEM time / kit time: above 1 the kit is faster)

## potential build (infinite)

| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |
|---|---|---|---|---|---|---|
| 64² | 0.626 (0.598-0.745) |  | 2.09 (1.76-2.7) | 1.68 (1.49-1.81) | 3.34x | 2.68x |
| 128² | 5.53 (3.13-7.33) |  | 9.69 (9.57-10.5) | 4.86 (4.23-7.11) | 1.75x | 0.88x |
| 192² | 13.2 (12.2-13.4) |  | 12.8 (12.2-13.5) | 9.97 (7.53-10.2) | 0.97x | 0.76x |
| 256² | 31.3 (21.1-35.8) |  | 41.4 (38.5-48.6) | 29.8 (27.9-32.1) | 1.32x | 0.95x |
| 512² | 155 (117-180) |  | 251 (220-379) | 120 (117-121) | 1.62x | 0.77x |

(ratios are abTEM time / kit time: above 1 the kit is faster)

## multislice 10 slices (plane wave, built potential)

| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |
|---|---|---|---|---|---|---|
| 64² | 0.929 (0.907-1.25) | 4.99 (4.49-5.39) | 3.49 (3.21-4.55) | 3.02 (2.76-3.12) | 3.76x | 3.25x |
| 128² | 11.6 (10.6-16.3) | 6.53 (5.5-10.5) | 17.5 (16.8-17.7) | 12 (11.1-12.6) | 1.50x | 1.84x |
| 192² | 25.9 (23.6-66.3) | 7.84 (6.74-8.13) | 20.4 (19.1-21.5) | 19.1 (18.7-20) | 0.79x | 2.44x |
| 256² | 60.1 (45.8-78.6) | 10.3 (9.34-10.7) | 71.4 (66.8-74.2) | 56.1 (53.7-65.6) | 1.19x | 5.42x |
| 512² | 299 (230-451) | 10.2 (10.1-11.1) | 678 (642-724) | 228 (225-237) | 2.27x | 22.43x |

(ratios are abTEM time / kit time: above 1 the kit is faster)

## transmission functions

| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |
|---|---|---|---|---|---|---|
| 64² | 0.527 (0.513-0.54) |  |  |  |  |  |
| 128² | 5.65 (4.63-15.6) |  |  |  |  |  |
| 192² | 13.6 (12-35.7) |  |  |  |  |  |
| 256² | 27.6 (23.4-65.1) |  |  |  |  |  |
| 512² | 107 (95.6-122) |  |  |  |  |  |

(ratios are abTEM time / kit time: above 1 the kit is faster)

## probe build

| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |
|---|---|---|---|---|---|---|
| 64² | 0.27 (0.161-0.296) |  | 0.84 (0.732-1.15) | 0.657 (0.587-1.09) | 3.11x | 2.43x |
| 128² | 0.824 (0.691-1.15) |  | 1.4 (1.28-1.66) | 1.62 (0.928-1.84) | 1.70x | 1.97x |
| 192² | 3.02 (2.09-4.85) |  | 2.19 (2.08-2.38) | 2.33 (1.92-2.72) | 0.73x | 0.77x |
| 256² | 3.64 (3.59-12.3) |  | 5.05 (4.32-5.49) | 3.04 (2.94-3.26) | 1.39x | 0.83x |
| 512² | 15.2 (14.4-56.8) |  | 42.4 (28.8-46.4) | 10.2 (9.46-10.8) | 2.78x | 0.67x |

(ratios are abTEM time / kit time: above 1 the kit is faster)

## STEM scan 8x8, annular

| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |
|---|---|---|---|---|---|---|
| 64² | 27.9 (19.5-39.8) | 9.39 (8.71-9.76) | 47.4 (47-51.7) | 33.8 (32-34.5) | 1.70x | 3.60x |
| 128² | 176 (143-246) | 13.8 (13.7-18.8) | 461 (451-464) | 206 (190-211) | 2.62x | 14.93x |
| 192² | 459 (386-611) | 35.3 (34.9-38.2) | 509 (493-559) | 367 (336-416) | 1.11x | 10.40x |
| 256² | 714 (448-768) | 50.6 (49.1-52.1) | 2.39e+03 (2.34e+03-2.64e+03) | 1.62e+03 (1.52e+03-1.66e+03) | 3.35x | 32.02x |
| 512² | 3.01e+03 (2.64e+03-3.34e+03) | 176 (174-179) | 8.84e+03 (8.75e+03-9.42e+03) | 8.39e+03 (8.1e+03-9.73e+03) | 2.94x | 47.55x |

(ratios are abTEM time / kit time: above 1 the kit is faster)

## frozen phonons 4 x multislice, built per configuration

| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |
|---|---|---|---|---|---|---|
| 64² | 10.6 (6.03-25.2) |  | 23.7 (20.9-24.6) | 21.5 (21.2-22.1) | 2.24x | 2.03x |
| 128² | 77.2 (58.9-119) |  | 130 (114-132) | 66.3 (61.2-78.9) | 1.68x | 0.86x |
| 192² | 173 (113-308) |  | 155 (138-220) | 94.9 (81.1-124) | 0.90x | 0.55x |
| 256² | 414 (314-559) |  | 472 (410-506) | 301 (290-458) | 1.14x | 0.73x |
| 512² | 1.32e+03 (1.12e+03-1.39e+03) |  | 1.64e+03 (1.62e+03-1.64e+03) | 1.38e+03 (1.34e+03-1.46e+03) | 1.24x | 1.04x |

(ratios are abTEM time / kit time: above 1 the kit is faster)

## multislice 10 slices (plane wave, potential built on the GPU)

| size | kit CPU f64 | kit Metal f32 | abTEM f64 | abTEM f32 | kit CPU / abTEM f64 | best kit / abTEM f32 |
|---|---|---|---|---|---|---|
| 64² |  | 6.04 (5.96-6.31) |  |  |  |  |
| 128² |  | 8.14 (7.63-8.89) |  |  |  |  |
| 192² |  | 10.1 (9.5-10.5) |  |  |  |  |
| 256² |  | 7.85 (6.79-11.7) |  |  |  |  |
| 512² |  | 15.7 (11.8-17.9) |  |  |  |  |

(ratios are abTEM time / kit time: above 1 the kit is faster)

## Backend decision per function

The smallest measured size at which Metal beats the CPU path. The break-even is bracketed by the sizes measured, not interpolated.

| function | decision | measured |
|---|---|---|
| fft2+ifft2 | CPU only: no Metal path yet (see open-items.md if it is a work item), or the work is CPU-side by nature | |
| potential build (infinite) | CPU only: no Metal path yet (see open-items.md if it is a work item), or the work is CPU-side by nature | |
| multislice 10 slices (plane wave, built potential) | Metal above 128² (break-even between 64² and 128²) | 64²: cpu 0.929 ms / metal 4.99 ms, 128²: cpu 11.6 ms / metal 6.53 ms, 192²: cpu 25.9 ms / metal 7.84 ms, 256²: cpu 60.1 ms / metal 10.3 ms, 512²: cpu 299 ms / metal 10.2 ms |
| transmission functions | CPU only: no Metal path yet (see open-items.md if it is a work item), or the work is CPU-side by nature | |
| probe build | CPU only: no Metal path yet (see open-items.md if it is a work item), or the work is CPU-side by nature | |
| STEM scan 8x8, annular | Metal above 64² (break-even between <64² and 64²) | 64²: cpu 27.9 ms / metal 9.39 ms, 128²: cpu 176 ms / metal 13.8 ms, 192²: cpu 459 ms / metal 35.3 ms, 256²: cpu 714 ms / metal 50.6 ms, 512²: cpu 3.01e+03 ms / metal 176 ms |
| frozen phonons 4 x multislice, built per configuration | CPU only: no Metal path yet (see open-items.md if it is a work item), or the work is CPU-side by nature | |
| multislice 10 slices (plane wave, potential built on the GPU) | Metal only (the GPU build of the Metal row above; its CPU counterpart is potential build + multislice) | |
