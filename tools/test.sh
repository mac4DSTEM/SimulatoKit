#!/bin/zsh
# SimulatoKit test gate: `swift test -c release`.
#
# Release, as SpectroKit's gate: the manifest carries no Debug `-O` (a versioned dependency cannot carry unsafe flags) and
# Debug at -Onone is many times slower on the FFT and potential loops; a consumer builds the kit optimised itself. Never run
# `swift build --build-tests -c release` (fails on this toolchain); this script is the gate.
#
# Stamp: a full run (no arguments) that passes, on a tree whose tracked Sources/, Tests/, Package.swift and Package.resolved
# are clean and hold no untracked file, writes .build/last-green (HEAD sha, UTC time). Any argument or a dirty tree writes nothing and says why on
# stderr. The stamp is local (.build/ is ignored); tools/last-green.sh prints it and whether HEAD still equals it.
set -euo pipefail
cd "${0:A:h}/.."
swift test -c release "$@"
if (( $# > 0 )); then
  print -u2 "tools/test.sh: green, but arguments were given (not a full run): no stamp"
elif [[ -n $(git status --porcelain --untracked-files=no -- Sources Tests Examples Package.swift Package.resolved) ]]; then
  print -u2 "tools/test.sh: green, but Sources/Tests/Examples/Package files are not clean: no stamp"
elif [[ -n $(git ls-files --others --exclude-standard -- Sources Tests Examples Package.swift Package.resolved) ]]; then
  # An untracked source or test file took part in the run but is not in HEAD (CuboKit's finding, 2026-10-15).
  print -u2 "tools/test.sh: green, but untracked files are in Sources/Tests/Examples: no stamp"
else
  mkdir -p .build
  print -r -- "commit $(git rev-parse HEAD)
time $(date -u +%Y-%m-%dT%H:%M:%SZ)" > .build/last-green
  print -u2 "tools/test.sh: stamped .build/last-green at $(git rev-parse --short HEAD)"
fi
