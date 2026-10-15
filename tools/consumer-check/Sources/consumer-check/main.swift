//
//  main.swift
//  Role: runs every section of consumer-check and exits non-zero on any failure. Fast: small grids, one second or two.
//

import Foundation

var report = Report()
checkCore(&report)
checkMultislice(&report)
checkProbe(&report)
checkDetector(&report)
exit(report.finish())
