# Consumers of SimulatoKit 0.2.1

Additive, no number moves: `MetalAccuracy` gains `detectorFloor` (1e−5) and `detectorBound(imageMaximum:incidentIntensity:)`,
the Metal detector bound accepted on 2026-10-15 (2e−4 × max(the image's maximum, 1e−5 × the incident intensity);
`docs/registration-performance.md`, amendment 4); a consumer that compares Metal detector images with CPU fixtures can use it.
