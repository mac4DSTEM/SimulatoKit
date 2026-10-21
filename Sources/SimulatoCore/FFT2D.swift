//
//  FFT2D.swift
//  Role: the CPU reference FFT: a complex 2-D DFT in float64 on Accelerate (vDSP_DFT_zop, double), NumPy's conventions:
//        the forward transform unnormalised with e^{−2πi k·r}, the inverse with 1/(nx ny). Rows (along y, contiguous)
//        then columns (along x, gathered). Any length: vDSP's DFT where it has one (f · 2^m, f in 1, 3, 5, 15, m ≥ 3),
//        Bluestein's chirp-z algorithm on a vDSP power-of-two convolution otherwise (own implementation, so the kit
//        stays standalone; ImagoKit's FFT is not a dependency).
//

import Accelerate

/// A 1-D complex DFT of one length in float64, unnormalised in both directions (the caller scales the inverse).
private final class LineDFT {
    let n: Int
    private let forward: vDSP_DFT_SetupD?
    private let inverse: vDSP_DFT_SetupD?
    // Bluestein: chirp w[j] = e^{−iπ j²/n}, the transform of conj(w) padded to m, and m-point work arrays (raw memory:
    // the loops below run per line, and bounds checks showed up in the profile).
    private let m: Int
    private let chirpRe: UnsafeMutablePointer<Double>
    private let chirpIm: UnsafeMutablePointer<Double>
    private let kernelRe: UnsafeMutablePointer<Double>
    private let kernelIm: UnsafeMutablePointer<Double>
    private let workRe: UnsafeMutablePointer<Double>
    private let workIm: UnsafeMutablePointer<Double>
    private let tmpRe: UnsafeMutablePointer<Double>
    private let tmpIm: UnsafeMutablePointer<Double>
    private let mForward: vDSP_DFT_SetupD?
    private let mInverse: vDSP_DFT_SetupD?

    init(n: Int) throws {
        guard n >= 1 else { throw SimulatoError.unsupportedGridLength(n) }
        self.n = n
        let direct = Grid.hasDirectTransform(n)
        var size = 16
        if !direct { while size < 2 * n - 1 { size *= 2 } }
        m = direct ? 0 : size
        let alloc = { (count: Int) -> UnsafeMutablePointer<Double> in
            let p = UnsafeMutablePointer<Double>.allocate(capacity: max(count, 1))
            p.initialize(repeating: 0, count: max(count, 1))
            return p
        }
        chirpRe = alloc(direct ? 0 : n)
        chirpIm = alloc(direct ? 0 : n)
        kernelRe = alloc(m)
        kernelIm = alloc(m)
        workRe = alloc(m)
        workIm = alloc(m)
        tmpRe = alloc(m)
        tmpIm = alloc(m)
        if direct {
            forward = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(n), .FORWARD)
            inverse = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(n), .INVERSE)
            mForward = nil
            mInverse = nil
            guard forward != nil, inverse != nil else { throw SimulatoError.unsupportedGridLength(n) }
            return
        }
        forward = nil
        inverse = nil
        mForward = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(size), .FORWARD)
        mInverse = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(size), .INVERSE)
        guard let mf = mForward, mInverse != nil else { throw SimulatoError.unsupportedGridLength(n) }
        for j in 0..<n {
            // j² mod 2n keeps the argument small and the chirp accurate for large n.
            let r = (j * j) % (2 * n)
            let angle = -Double.pi * Double(r) / Double(n)
            chirpRe[j] = cos(angle)
            chirpIm[j] = sin(angle)
            kernelRe[j] = chirpRe[j]
            kernelIm[j] = -chirpIm[j]
            if j > 0 {
                kernelRe[size - j] = chirpRe[j]
                kernelIm[size - j] = -chirpIm[j]
            }
        }
        vDSP_DFT_ExecuteD(mf, kernelRe, kernelIm, tmpRe, tmpIm)
        kernelRe.update(from: tmpRe, count: size)
        kernelIm.update(from: tmpIm, count: size)
    }

    deinit {
        if let forward { vDSP_DFT_DestroySetupD(forward) }
        if let inverse { vDSP_DFT_DestroySetupD(inverse) }
        if let mForward { vDSP_DFT_DestroySetupD(mForward) }
        if let mInverse { vDSP_DFT_DestroySetupD(mInverse) }
        for p in [chirpRe, chirpIm, kernelRe, kernelIm, workRe, workIm, tmpRe, tmpIm] { p.deallocate() }
    }

    /// `re`/`im` hold n values and are overwritten by their DFT (e^{−2πi jk/n}, or e^{+…} when `inverse`), unscaled.
    /// The scratch pointers hold at least n values (used by the direct transform only).
    func execute(_ re: UnsafeMutablePointer<Double>, _ im: UnsafeMutablePointer<Double>, inverse isInverse: Bool,
                 scratchRe: UnsafeMutablePointer<Double>, scratchIm: UnsafeMutablePointer<Double>) {
        if let forward, let inverse {
            vDSP_DFT_ExecuteD(isInverse ? inverse : forward, re, im, scratchRe, scratchIm)
            re.update(from: scratchRe, count: n)
            im.update(from: scratchIm, count: n)
            return
        }
        // The inverse is conj(DFT(conj x)).
        let sign = isInverse ? -1.0 : 1.0
        for j in 0..<n {
            let xr = re[j], xi = sign * im[j]
            workRe[j] = xr * chirpRe[j] - xi * chirpIm[j]
            workIm[j] = xr * chirpIm[j] + xi * chirpRe[j]
        }
        (workRe + n).update(repeating: 0, count: m - n)
        (workIm + n).update(repeating: 0, count: m - n)
        vDSP_DFT_ExecuteD(mForward!, workRe, workIm, tmpRe, tmpIm)
        for j in 0..<m {
            let a = tmpRe[j], b = tmpIm[j]
            workRe[j] = a * kernelRe[j] - b * kernelIm[j]
            workIm[j] = a * kernelIm[j] + b * kernelRe[j]
        }
        vDSP_DFT_ExecuteD(mInverse!, workRe, workIm, tmpRe, tmpIm)
        let scale = 1.0 / Double(m)
        for k in 0..<n {
            let a = tmpRe[k] * scale, b = tmpIm[k] * scale
            re[k] = a * chirpRe[k] - b * chirpIm[k]
            im[k] = sign * (a * chirpIm[k] + b * chirpRe[k])
        }
    }
}

/// A 2-D complex DFT plan for one grid size. Not shareable across threads; make one per task.
package final class FFT2D {
    package let nx: Int
    package let ny: Int
    private let planX: LineDFT
    private let planY: LineDFT
    private var lineRe: [Double]
    private var lineIm: [Double]
    private var outRe: [Double]
    private var outIm: [Double]

    package init(nx: Int, ny: Int) throws {
        guard Grid.isSupportedLength(nx) else { throw SimulatoError.unsupportedGridLength(nx) }
        guard Grid.isSupportedLength(ny) else { throw SimulatoError.unsupportedGridLength(ny) }
        self.nx = nx
        self.ny = ny
        planX = try LineDFT(n: nx)
        planY = nx == ny ? planX : try LineDFT(n: ny)
        let n = max(nx, ny)
        lineRe = [Double](repeating: 0, count: n)
        lineIm = [Double](repeating: 0, count: n)
        outRe = [Double](repeating: 0, count: n)
        outIm = [Double](repeating: 0, count: n)
    }

    /// In place, unnormalised, e^{−2πi k·r} (NumPy's `fft2`).
    package func forward(_ field: inout ComplexField) {
        transform(&field.real, &field.imaginary, inverse: false)
    }

    /// In place, with 1/(nx ny) (NumPy's `ifft2`).
    package func inverse(_ field: inout ComplexField) {
        transform(&field.real, &field.imaginary, inverse: true)
        let scale = 1.0 / Double(nx * ny)
        for i in field.real.indices {
            field.real[i] *= scale
            field.imaginary[i] *= scale
        }
    }

    private func transform(_ re: inout [Double], _ im: inout [Double], inverse: Bool) {
        precondition(re.count == nx * ny && im.count == nx * ny)
        let nx = self.nx, ny = self.ny
        re.withUnsafeMutableBufferPointer { r in
            im.withUnsafeMutableBufferPointer { m in
                outRe.withUnsafeMutableBufferPointer { orr in
                    outIm.withUnsafeMutableBufferPointer { oi in
                        // Along y: each row ix is contiguous.
                        for ix in 0..<nx {
                            let base = ix * ny
                            planY.execute(r.baseAddress! + base, m.baseAddress! + base, inverse: inverse,
                                          scratchRe: orr.baseAddress!, scratchIm: oi.baseAddress!)
                        }
                        lineRe.withUnsafeMutableBufferPointer { lr in
                            lineIm.withUnsafeMutableBufferPointer { li in
                                // Along x: gather column iy.
                                for iy in 0..<ny {
                                    for ix in 0..<nx {
                                        lr[ix] = r[ix * ny + iy]
                                        li[ix] = m[ix * ny + iy]
                                    }
                                    planX.execute(lr.baseAddress!, li.baseAddress!, inverse: inverse,
                                                  scratchRe: orr.baseAddress!, scratchIm: oi.baseAddress!)
                                    for ix in 0..<nx {
                                        r[ix * ny + iy] = lr[ix]
                                        m[ix * ny + iy] = li[ix]
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
