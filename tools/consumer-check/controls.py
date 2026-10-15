#!/usr/bin/env python3
"""consumer-check negative controls (run through `controls.sh`; Python 3 standard library only).

The consumer package (tools/consumer-check) compiles only while every symbol it names is `public`. This proves it: for each
of a few representative declarations across the four modules, ONE `public` is reverted to `package` in a COPY of the kit
(rsync into <root>/.build/consumer-check-controls/<kit name>, never the working tree), the consumer is built against the
copy (its Package.swift takes the dependency by path), and the build must FAIL with a compiler error that names that symbol.
A control is caught only when
  * the unedited copy builds first (the baseline, so a failure is the edit's), and
  * the edit applied exactly once (the old text occurs once in its file), and
  * the build exits non-zero and an `error:` line matches the control's pattern, which names the symbol:
      "access":  the consumer is refused, `'<member>' is inaccessible due to 'package' protection level`, or `cannot find
                 '<Type>' in scope` for a whole type;
      "overload": the reverted declaration is one of two overloads (an initialiser), so the compiler resolves the call to the
                 public sibling and says so, naming the reverted labels: `incorrect argument label in call (have 'a:b:', expected
                 'c:d:')`, or `extra arguments at positions ...` at a consumer line that holds the reverted labels (`token`);
      "library": the kit itself refuses the edit, because another public declaration is typed with the now-package one
                 (`... cannot be declared public because its parameter uses a package type`, at a source line that names the type).
One copy and one warm .build serve every control (the edited file is restored before the next), so each rebuild is incremental.

usage: controls.py <kit root> [<control name> ...]     (no name: all of them)
exit:  0 only when the baseline builds and every control is caught.
"""
import os
import re
import subprocess
import sys
import time

ROOT = os.path.abspath(sys.argv[1])
ONLY = sys.argv[2:]
NAME = os.path.basename(ROOT)
BASE = os.path.join(ROOT, ".build", "consumer-check-controls")
COPY = os.path.join(BASE, NAME)
CONSUMER = os.path.join(COPY, "tools", "consumer-check")
ANSI = re.compile(r"\x1b\[[0-9;]*m")

MS = "Sources/SimulatoMultislice/"
PR = "Sources/SimulatoProbe/"
DT = "Sources/SimulatoDetector/"
CO = "Sources/SimulatoCore/"


def control(name, rel, old, pattern, token=None):
    """`old` (a `public ...` declaration, once in `rel`) becomes `package ...`; the build must fail with `pattern` in an error line
    (and, with `token`, at a source line that contains `token`, for a diagnostic that does not itself name the symbol)."""
    assert old.startswith("public ") or " public " in old, name
    return dict(name=name, rel=rel, old=old, new=old.replace("public ", "package ", 1), pattern=pattern, token=token)


def inaccessible(member):
    return r"'%s' is inaccessible due to 'package' protection level" % re.escape(member)


CONTROLS = [
    control("ProjectedPotential.init(structure:nx:ny:sliceThickness:parameterization:)", MS + "ProjectedPotential.swift",
            "public init(structure: AtomicStructure, nx: Int, ny: Int, sliceThickness: Double,",
            r"incorrect argument label in call \(have 'structure:nx:ny:sliceThickness:'"),
    control("Multislice.run", MS + "Multislice.swift",
            "public static func run(_ incident: Waves, through potential: ProjectedPotential",
            # With the ProjectedPotential overload hidden, the Potential overload is chosen and the argument is refused.
            inaccessible("run") + r"|cannot convert value of type 'ProjectedPotential' to expected argument type 'Potential'"),
    control("Potential", MS + "Potential.swift", "public struct Potential: Sendable, Hashable {",
            r"cannot find 'Potential' in scope|'Potential' is inaccessible due to 'package' protection level"
            r"|cannot be declared public because its parameter uses a package type"),
    control("Probe.build(at:)", PR + "Probe.swift", "public func build(at position: SIMD2<Double>)",
            inaccessible("build")),
    control("AnnularDetector", DT + "Detectors.swift", "public struct AnnularDetector: Detector, Hashable {",
            r"cannot find 'AnnularDetector' in scope|'AnnularDetector' is inaccessible due to 'package' protection level"),
    control("ContrastTransferFunction", PR + "ContrastTransferFunction.swift",
            "public struct ContrastTransferFunction: Sendable, Hashable {",
            r"cannot be declared public because .* uses a package type", token="ContrastTransferFunction"),
    control("GridScan.positions", PR + "GridScan.swift", "public var positions: [SIMD2<Double>] {",
            inaccessible("positions")),
    control("FrozenPhonons.init(structure:sigmas:configurations:seed:directions:)", MS + "FrozenPhonons.swift",
            "public init(structure: AtomicStructure, sigmas: [Int: Double], configurations: Int, seed: UInt64,",
            r"extra arguments at positions #2, #3, #4 in call", token="sigmas:"),
    control("Grid.fftfreq", CO + "Grid.swift", "public static func fftfreq(_ n: Int, _ d: Double) -> [Double] {",
            inaccessible("fftfreq")),
    control("Detector (protocol)", DT + "Detectors.swift", "public protocol Detector: Sendable {",
            r"cannot be declared public because .* uses a package type", token="Detector"),
]


def run(cmd, cwd):
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return p.returncode, ANSI.sub("", p.stdout)


def sync():
    """The kit into the copy, without .build, References, .git; the copy's own .build (the warm one) stays."""
    os.makedirs(COPY, exist_ok=True)
    code, out = run(["rsync", "-a", "--delete", "--exclude", ".build", "--exclude", "References", "--exclude", ".git",
                     ROOT + "/", COPY + "/"], ROOT)
    if code != 0:
        sys.exit("rsync failed:\n" + out)


def build():
    return run(["swift", "build", "-c", "release", "-j", "4"], CONSUMER)


def error_site(line):
    """The source text of the line an error is reported at (`path:line:col: error: ...`), or ''."""
    m = re.match(r"^(/\S+?\.swift):(\d+):\d+: error:", line.strip())
    if not m:
        return ""
    try:
        with open(m.group(1)) as f:
            return f.read().splitlines()[int(m.group(2)) - 1]
    except (OSError, IndexError):
        return ""


def first_errors(out, n=2):
    return [l.strip() for l in out.splitlines() if "error:" in l][:n]


def main():
    started = time.time()
    sync()
    t0 = time.time()
    code, out = build()
    ok = code == 0
    print("%s  baseline: the unedited copy builds  (%.0f s)" % ("CAUGHT" if ok else "MISSED", time.time() - t0))
    if not ok:
        print("\n".join(first_errors(out, 5)))
        sys.exit(1)
    missed = 0
    selected = [c for c in CONTROLS if not ONLY or c["name"] in ONLY]
    for c in selected:
        path = os.path.join(COPY, c["rel"])
        with open(path) as f:
            original = f.read()
        if original.count(c["old"]) != 1:
            print("MISSED  %s: the declaration text occurs %d times in %s (want 1)" % (c["name"], original.count(c["old"]), c["rel"]))
            missed += 1
            continue
        t0 = time.time()
        try:
            with open(path, "w") as f:
                f.write(original.replace(c["old"], c["new"], 1))
            code, out = build()
        finally:
            with open(path, "w") as f:
                f.write(original)
        errors = [l for l in out.splitlines() if "error:" in l]
        matching = [l for l in errors if re.search(c["pattern"], l) and (not c["token"] or c["token"] in error_site(l))]
        hit = code != 0 and bool(matching)
        if hit:
            line = matching[0]
            print("CAUGHT  %s  (%.0f s)\n          %s" % (c["name"], time.time() - t0, re.sub(r"^.*/Sources/", "Sources/", line.strip())))
        else:
            missed += 1
            why = "the build succeeded" if code == 0 else "no error line matched /%s/%s" % (c["pattern"], " at a line with " + c["token"] if c["token"] else "")
            print("MISSED  %s: %s" % (c["name"], why))
            print("\n".join("          " + l for l in first_errors(out, 3)))
    # the copy is restored: the unedited copy must still build (the restore did not leave an edit behind)
    t0 = time.time()
    code, out = build()
    print("%s  after the controls: the restored copy builds again  (%.0f s)" % ("CAUGHT" if code == 0 else "MISSED", time.time() - t0))
    if code != 0:
        missed += 1
    print("consumer-check controls: %d of %d caught (+ baseline and restore), %.0f s" % (len(selected) - missed, len(selected), time.time() - started))
    sys.exit(0 if missed == 0 else 1)


main()
