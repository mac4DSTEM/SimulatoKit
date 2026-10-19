#!/bin/zsh
# SimulatoKit README consumer: proves the README's code alone is enough for a stranger. Extracts every ```swift block
# from the "## Quick start" heading on, in order, VERBATIM into readme/Sources/readme-consumer/main.swift,
# appends readme/checks.swift.txt (assertions on the README's own `hrtem`, `haadf` and `scan`), builds the separate package
# readme/ in RELEASE (-j 4, no warnings allowed) and runs it. Exit 0 only if all of that holds.
#
#   readme.sh   prints "readme-consumer: the README's quick start builds and runs (N checks)" and exits 0.
set -euo pipefail
here=${0:A:h}
root=${here:h:h}
pkg=$here/readme
die() { print -u2 "readme-consumer: $*"; exit 1 }

main=$pkg/Sources/readme-consumer/main.swift
mkdir -p ${main:h}
# every ```swift block from the "## Quick start" heading on, without the fences, in order
awk '/^## Quick start[[:space:]]*$/ {s=1; next} s && !f && /^```swift[[:space:]]*$/ {f=1; next} f && /^```[[:space:]]*$/ {f=0; next} f {print}' \
  $root/README.md >$main.code
[[ -s $main.code ]] || die "no \`\`\`swift block found after the '## Quick start' heading of README.md"
for name in hrtem haadf scan; do
  grep -qE "^(let|var) $name\\b" $main.code \
    || die "the README's quick start defines no top-level 'let/var $name' (the checks need hrtem, haadf and scan)"
done
# top-level code in main.swift may `try` (Swift 6), so the README code needs no wrapper
cat $main.code $pkg/checks.swift.txt >$main
rm $main.code

cd $pkg
mkdir -p .build
log=.build/readme-consumer-build.log
# --disable-dependency-cache: SwiftPM's shared repository cache can print a 'skipping cache' warning on a CI runner,
# which the warning count would take for the kit's own (DiffractoKit 0.4.1's CI).
if ! swift build --disable-dependency-cache -c release -j 4 >$log 2>&1; then
  grep -E "error" $log | head -20 >&2
  die "BUILD FAILED (full log: $pkg/$log)"
fi
warnings=$(grep -c "warning:" $log || true)
(( warnings == 0 )) || { grep "warning:" $log | head >&2; die "$warnings warning(s) in the build" }
bin=$(swift build --disable-dependency-cache -c release --show-bin-path)
out=.build/readme-consumer-run.log
set +e
$bin/readme-consumer >$out 2>&1
rc=$?
set -e
cat $out
(( rc == 0 )) || die "FAILED (exit $rc)"
summary=$(grep "^readme-consumer summary:" $out) || die "no summary line in the run output"
checks=${${summary#readme-consumer summary: }%% checks*}
print "readme-consumer: the README's quick start builds and runs ($checks checks)"
