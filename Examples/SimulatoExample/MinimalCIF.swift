//
//  MinimalCIF.swift
//  Role: the example's CIF reader: cell, symmetry operations and atom sites of an orthogonal cell, expanded into
//        the cell's atoms. It is not part of the library (the kit stays standalone,
//        docs/adr/0001-standalone-structure-types.md) and makes no claim beyond the example (no uncertainty handling
//        beyond stripping "(n)", no written-precision merging, non-orthogonal cells refused).
//

import Foundation
import SimulatoCore

enum MinimalCIF {
    struct Failure: Error, CustomStringConvertible {
        let description: String
    }

    static let symbols = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
                          "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
                          "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
                          "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
                          "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
                          "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
                          "Es", "Fm", "Md"]

    /// The unit cell's atoms in an orthogonal box.
    static func read(_ text: String) throws -> (name: String, cell: AtomicStructure) {
        let tokens = tokenize(text)
        var scalars: [String: String] = [:]
        var loops: [[String: [String]]] = []
        var name = "structure"
        var i = 0
        while i < tokens.count {
            let t = tokens[i]
            if t.lowercased().hasPrefix("data_") { name = String(t.dropFirst(5)); i += 1; continue }
            if t.lowercased() == "loop_" {
                var tags: [String] = []
                i += 1
                while i < tokens.count, tokens[i].hasPrefix("_") { tags.append(tokens[i].lowercased()); i += 1 }
                var values: [String] = []
                while i < tokens.count, !tokens[i].hasPrefix("_"), tokens[i].lowercased() != "loop_",
                      !tokens[i].lowercased().hasPrefix("data_") {
                    values.append(tokens[i]); i += 1
                }
                guard !tags.isEmpty, values.count % tags.count == 0 else { throw Failure(description: "truncated loop") }
                var columns: [String: [String]] = [:]
                for (k, tag) in tags.enumerated() {
                    columns[tag] = stride(from: k, to: values.count, by: tags.count).map { values[$0] }
                }
                loops.append(columns)
                continue
            }
            if t.hasPrefix("_"), i + 1 < tokens.count { scalars[t.lowercased()] = tokens[i + 1]; i += 2; continue }
            i += 1
        }
        func number(_ tag: String) throws -> Double {
            guard let s = scalars[tag], let v = Double(s.split(separator: "(")[0]) else {
                throw Failure(description: "missing or unparsable \(tag)")
            }
            return v
        }
        let a = try number("_cell_length_a"), b = try number("_cell_length_b"), c = try number("_cell_length_c")
        for tag in ["_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"] where abs(try number(tag) - 90) > 1e-6 {
            throw Failure(description: "\(tag) is not 90°: the example simulates orthogonal cells only")
        }
        let opTags = ["_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"]
        let ops = loops.compactMap { loop in opTags.compactMap { loop[$0] }.first }.first ?? ["x,y,z"]
        guard let sites = loops.first(where: { $0["_atom_site_fract_x"] != nil }) else {
            throw Failure(description: "no _atom_site_fract_x loop")
        }
        let labels = sites["_atom_site_type_symbol"] ?? sites["_atom_site_label"] ?? []
        func frac(_ s: String) -> Double? { Double(s.split(separator: "(")[0]) }
        var atoms: [AtomicStructure.Atom] = []
        var seen: [(Int, SIMD3<Double>)] = []
        for k in labels.indices {
            let letters = labels[k].prefix { $0.isLetter }
            let symbol = letters.count >= 2 && symbols.contains(String(letters.prefix(2))) ? String(letters.prefix(2))
                : String(letters.prefix(1))
            guard let z = symbols.firstIndex(of: symbol).map({ $0 + 1 }) else {
                throw Failure(description: "unknown element in \(labels[k])")
            }
            guard let x = frac(sites["_atom_site_fract_x"]![k]), let y = frac(sites["_atom_site_fract_y"]![k]),
                  let w = frac(sites["_atom_site_fract_z"]![k]) else { throw Failure(description: "unparsable site \(k)") }
            for op in ops {
                var p = try apply(op, SIMD3(x, y, w))
                for d in 0..<3 { p[d] -= p[d].rounded(.down); if p[d] > 1 - 1e-9 { p[d] = 0 } }
                if seen.contains(where: { $0.0 == z && distance($0.1, p) < 1e-3 }) { continue }
                seen.append((z, p))
                atoms.append(.init(atomicNumber: z, position: SIMD3(p.x * a, p.y * b, p.z * c)))
            }
        }
        return (name, try AtomicStructure(box: SIMD3(a, b, c), atoms: atoms))
    }

    /// Distance in fractional coordinates, periodic.
    static func distance(_ p: SIMD3<Double>, _ q: SIMD3<Double>) -> Double {
        var d = p - q
        for k in 0..<3 { d[k] -= d[k].rounded() }
        return (d * d).sum().squareRoot()
    }

    /// One `x,y,z`-form operation applied to a fractional position.
    static func apply(_ op: String, _ p: SIMD3<Double>) throws -> SIMD3<Double> {
        let parts = op.replacingOccurrences(of: " ", with: "").lowercased().split(separator: ",")
        guard parts.count == 3 else { throw Failure(description: "unparsable operation \(op)") }
        var out = SIMD3<Double>()
        for (k, part) in parts.enumerated() {
            var value = 0.0
            var term = ""
            for ch in part + "+" {
                if (ch == "+" || ch == "-") && !term.isEmpty && term != "+" && term != "-" {
                    value += try evaluate(term, p)
                    term = String(ch)
                } else {
                    term.append(ch)
                }
            }
            out[k] = value
        }
        return out
    }

    static func evaluate(_ term: String, _ p: SIMD3<Double>) throws -> Double {
        var sign = 1.0, t = Substring(term)
        if t.first == "-" { sign = -1; t = t.dropFirst() } else if t.first == "+" { t = t.dropFirst() }
        for (axis, letter) in ["x", "y", "z"].enumerated() where t.hasSuffix(letter) {
            let coefficient = t.dropLast()
            let c = coefficient.isEmpty ? 1.0 : try fraction(coefficient.hasSuffix("*") ? coefficient.dropLast() : coefficient)
            return sign * c * p[axis]
        }
        return sign * (try fraction(t))
    }

    static func fraction(_ s: Substring) throws -> Double {
        let f = s.split(separator: "/")
        if f.count == 2, let n = Double(f[0]), let d = Double(f[1]) { return n / d }
        if f.count == 1, let v = Double(f[0]) { return v }
        throw Failure(description: "unparsable number \(s)")
    }

    /// CIF tokens: whitespace-separated, quotes kept together, comments dropped, semicolon text fields skipped.
    static func tokenize(_ text: String) -> [String] {
        var tokens: [String] = []
        var inText = false
        for raw in text.components(separatedBy: .newlines) {
            if raw.hasPrefix(";") { inText.toggle(); continue }
            if inText { continue }
            var current = "", quote: Character? = nil
            for ch in raw {
                if let q = quote {
                    if ch == q { quote = nil; tokens.append(current); current = "" } else { current.append(ch) }
                } else if ch == "#" {
                    break
                } else if ch == "'" || ch == "\"" {
                    quote = ch
                } else if ch == " " || ch == "\t" {
                    if !current.isEmpty { tokens.append(current); current = "" }
                } else {
                    current.append(ch)
                }
            }
            if !current.isEmpty { tokens.append(current) }
        }
        return tokens
    }
}
