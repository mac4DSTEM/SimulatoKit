//
//  Report.swift
//  Role: the check ledger of consumer-check: one printed line per check, a count of failures, no test framework (a
//        consumer package, not a test target). `finish()` is the process exit status.
//

import Foundation
import SimulatoCore

struct Report {
    private(set) var passed = 0
    private(set) var failed = 0
    private(set) var skipped = 0

    mutating func check(_ name: String, _ ok: Bool, _ detail: @autoclosure () -> String = "") {
        if ok { passed += 1; print("PASS  \(name)") } else { failed += 1; print("FAIL  \(name)  [\(detail())]") }
    }
    mutating func equal<T: Equatable>(_ name: String, _ got: T, _ want: T) {
        check(name, got == want, "got \(got), want \(want)")
    }
    mutating func near(_ name: String, _ got: Double, _ want: Double, _ tolerance: Double) {
        check(name, abs(got - want) <= tolerance, "got \(got), want \(want) +- \(tolerance)")
    }
    mutating func nearRelative(_ name: String, _ got: Double, _ want: Double, _ relative: Double) {
        check(name, abs(got - want) <= relative * abs(want), "got \(got), want \(want) rel \(relative)")
    }
    /// A named skip: printed with its reason, not counted as a check.
    mutating func skip(_ name: String, _ reason: String) {
        skipped += 1
        print("SKIP  \(name): \(reason)")
    }
    /// The body must throw exactly `expected` (a refusal of the kit).
    mutating func refuses<T>(_ name: String, _ expected: SimulatoError, _ body: () throws -> T) {
        do {
            _ = try body()
            check(name, false, "no error thrown, wanted \(expected)")
        } catch {
            check(name, (error as? SimulatoError) == expected, "threw \(error), wanted \(expected)")
        }
    }
    /// The body must throw a `SimulatoError` of the named case (payload not compared).
    mutating func refusesCase<T>(_ name: String, _ caseName: String, _ body: () throws -> T) {
        do {
            _ = try body()
            check(name, false, "no error thrown, wanted \(caseName)")
        } catch {
            let label = Mirror(reflecting: error).children.first?.label ?? "\(error)"
            check(name, error is SimulatoError && label == caseName, "threw \(error) (\(label)), wanted \(caseName)")
        }
    }
    /// Runs `body`; a thrown error is one failed check, not a crash, so every later check still runs.
    mutating func section(_ name: String, _ body: (inout Report) throws -> Void) {
        var copy = self
        do { try body(&copy) } catch { copy.check("\(name): ran to the end", false, "threw \(error)") }
        self = copy
    }
    func finish() -> Int32 {
        print("consumer-check summary: \(passed + failed) checks, \(passed) passed, \(failed) failed, \(skipped) skipped")
        return failed == 0 ? 0 : 1
    }
}

/// The value, or a thrown check failure (a consumer package has no #require).
func unwrap<T>(_ value: T?, _ what: String = #function) throws -> T {
    guard let value else { throw CheckFailure("nil where a value was expected (\(what))") }
    return value
}

struct CheckFailure: Error, CustomStringConvertible {
    let description: String
    init(_ description: String) { self.description = description }
}

/// The largest |a − b| over two arrays of equal length.
func maxDifference(_ a: [Double], _ b: [Double]) -> Double {
    precondition(a.count == b.count)
    var m = 0.0
    for i in a.indices { m = max(m, abs(a[i] - b[i])) }
    return m
}
