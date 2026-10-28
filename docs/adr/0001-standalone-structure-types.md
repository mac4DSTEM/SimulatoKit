# 1. SimulatoKit keeps its own structure types

Date: 2026-10-28

## Status

Accepted.

## Context

SimulatoCore needs atoms in a box to slice a potential. DiffractoKit, a sibling library, has richer crystal types
(`CrystalStructure`, a CIF reader). Earlier comments marked `AtomicStructure` and the example's `MinimalCIF` "for
replacement" by DiffractoKit's types once DiffractoKit was tagged. DiffractoKit is tagged now, and the kits are kept
standalone: each depends on no other kit, so each can be used, versioned and released on its own.

## Decision

`AtomicStructure` (an orthogonal box with Cartesian atoms) stays in SimulatoCore, and `MinimalCIF` stays in the example.
SimulatoKit takes no dependency on DiffractoKit. An app that holds a DiffractoKit `CrystalStructure` converts it to an
`AtomicStructure` itself.

## Consequences

- Non-orthogonal and rotated cells, `orthogonalize_cell` and ASE's `pbc` are SimulatoKit's own work items, not a
  dependency to wait for (`docs/parity.md`).
- `MinimalCIF` stays a small, example-only reader; it is not part of the library.
- Removing either type, or adding a dependency on another kit, reverses this decision and needs a new ADR.
