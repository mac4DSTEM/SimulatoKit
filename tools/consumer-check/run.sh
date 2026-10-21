#!/bin/zsh
# SimulatoKit consumer check: build tools/consumer-check (a separate package that depends on SimulatoKit by path) in
# RELEASE, as the kit's tests run (the kit's tests use -c release), and run it. Exit 0 only if the build is
# warning-free, no test-only target was built, and every check passed.
#
#   run.sh            build, run the checks, print "consumer-check: all passed — N checks ..." and exit 0.
#
set -euo pipefail
here=${0:A:h}
root=${here:h:h}
die() { print -u2 "consumer-check: $*"; exit 1 }
cd $here
mkdir -p .build
log=.build/consumer-check-build.log
# --disable-dependency-cache: SwiftPM's shared repository cache can print a 'skipping cache' warning on a CI runner,
# which the warning count would take for the kit's own (DiffractoKit 0.4.1's CI).
if ! swift build --disable-dependency-cache -c release -j 4 >$log 2>&1; then
  grep -E "error" $log | head -20 >&2
  die "BUILD FAILED (full log: $here/$log)"
fi
warnings=$(grep -c "warning:" $log || true)
(( warnings == 0 )) || { grep "warning:" $log | head >&2; die "$warnings warning(s) in the build" }
if grep -qE "SimulatoCoreTests|SimulatoMultisliceTests|SimulatoProbeTests|SimulatoDetectorTests|SimulatoKitParityTests|SimulatoKitTestSupport|SimulatoExample" $log; then
  die "the build log names a test-only or example target"
fi
bin=$(swift build --disable-dependency-cache -c release --show-bin-path)
out=.build/consumer-check-run.log
set +e
$bin/consumer-check >$out 2>&1
rc=$?
set -e
cat $out
(( rc == 0 )) || die "FAILED (exit $rc)"
summary=$(grep "^consumer-check summary:" $out)
checks=${${summary#consumer-check summary: }%% checks*}
head=$(git -C $root --no-optional-locks rev-parse --short HEAD 2>/dev/null || print unknown)
print "consumer-check: all passed — $checks checks through SimulatoKit's public API only, release build, no warnings, no test-only target built (SimulatoKit $head)"
# the second consumer: the README's Quick start, extracted verbatim and run (readme.sh)
$here/readme.sh
