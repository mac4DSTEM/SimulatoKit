#!/bin/zsh
# SimulatoKit consumer-check negative controls: reverts one `public` to `package` in a copy of the kit, nine times, and
# requires the consumer to fail to build with an error naming the symbol (see controls.py). A cold copy build, then
# incremental ones: about a minute. Exit 0 only if the baseline builds and every control is caught.
#
#   controls.sh                 all controls
#   controls.sh "Probe.build(at:)" ...    only the named ones
set -euo pipefail
here=${0:A:h}
exec /usr/bin/python3 -I $here/controls.py ${here:h:h} "$@"
