#!/bin/bash
# Proof for no-sleep-poll.sh: feeds synthetic hook JSON, prints PASS/FAIL per case, exits non-zero on any failure.
cd "$(dirname "$0")" || exit 1
fail=0
check() { # expected-exit command
  local json; json=$(/usr/bin/python3 -I -c 'import json,sys; print(json.dumps({"tool_name":"Bash","tool_input":{"command":sys.argv[1]}}))' "$2")
  ./no-sleep-poll.sh <<<"$json" 2>/dev/null; local rc=$?
  if [[ $rc == "$1" ]]; then echo "PASS (exit $rc): $2"; else echo "FAIL (exit $rc, wanted $1): $2"; fail=1; fi
}
check 2 'sleep 60'
check 2 'sleep 30'
check 2 'sleep 45 && echo done'
check 2 'echo a; /bin/sleep 120'
check 2 'sleep 1m'
check 2 'sleep 0.5h'
check 2 'while true; do x; sleep 2; done'
check 2 'until grep -q x f; do sleep 1; done'
check 2 'pkill -f swift-test'
check 2 'pkill -9f foo'
check 2 'pkill --full x'
check 2 'killall swift'
check 2 'cd x && killall -9 y'
check 2 'echo a; /usr/bin/pkill -f x'
check 2 'git log --grep killall'   # matched as text: refused (accepted false positive)
check 0 'echo "killall swift"'   # quoted right after a quote: passes (same text policy as sleep)
check 0 'kill 4242'
check 0 'kill -9 4242'
check 0 'pkill -P 4242'
check 0 'pkill swift-test'   # no -f: name match only, allowed
check 0 'pkill -9 -P 4242 && ls -f'   # -f belongs to ls, not pkill
check 0 'echo "sleep 60"'   # quoted: passes (sleep must follow whitespace or an operator)
check 0 'sleep 1'
check 0 'sleep 5 && ls'
check 0 'sleep 29.9'
check 0 'git log'
check 0 'swift test'
check 0 'tools/test.sh --filter Sleeper'
exit $fail
