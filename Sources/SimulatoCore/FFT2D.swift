//
//  FFT2D.swift
//  Role: the CPU reference FFT: a complex 2-D DFT in float64 on Accelerate (vDSP_DFT_zop, double), NumPy's conventions:
//        the forward transform unnormalised with e^{−2πi k·r}, the inverse with 1/(nx ny). Rows (along y, contiguous)
//        then columns (along x, gathered).
//

import Accelerate

/// A 2-D complex DFT plan for one grid size. Not shareable across threads; make one per task.
package final class FFT2D {
    package let nx: Int
    package let ny: Int
    private let forwardX: vDSP_DFT_SetupD
    private let inverseX: vDSP_DFT_SetupD
    private let forwardY: vDSP_DFT_SetupD
    private let inverseY: vDSP_DFT_SetupD
    private var lineRe: [Double]
    private var lineIm: [Double]
    private var outRe: [Double]
    private var outIm: [Double]

    package init(nx: Int, ny: Int) throws {
        guard Grid.isSupportedLength(nx) else { throw SimulatoError.unsupportedGridLength(nx) }
        guard Grid.isSupportedLength(ny) else { throw SimulatoError.unsupportedGridLength(ny) }
        guard let fx = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(nx), .FORWARD),
              let ix = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(nx), .INVERSE),
              let fy = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(ny), .FORWARD),
              let iy = vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(ny), .INVERSE)
        else { throw SimulatoError.unsupportedGridLength(nx) }
        self.nx = nx
        self.ny = ny
        forwardX = fx
        inverseX = ix
        forwardY = fy
        inverseY = iy
        let n = max(nx, ny)
        lineRe = [Double](repeating: 0, count: n)
        lineIm = [Double](repeating: 0, count: n)
        outRe = [Double](repeating: 0, count: n)
        outIm = [Double](repeating: 0, count: n)
    }

    deinit {
        vDSP_DFT_DestroySetupD(forwardX)
        vDSP_DFT_DestroySetupD(inverseX)
        vDSP_DFT_DestroySetupD(forwardY)
        vDSP_DFT_DestroySetupD(inverseY)
    }

    /// In place, unnormalised, e^{−2πi k·r} (NumPy's `fft2`).
    package func forward(_ field: inout ComplexField) {
        transform(&field.real, &field.imaginary, rows: forwardY, columns: forwardX)
    }

    /// In place, with 1/(nx ny) (NumPy's `ifft2`).
    package func inverse(_ field: inout ComplexField) {
        transform(&field.real, &field.imaginary, rows: inverseY, columns: inverseX)
        let scale = 1.0 / Double(nx * ny)
        for i in field.real.indices {
            field.real[i] *= scale
            field.imaginary[i] *= scale
        }
    }

    private func transform(_ re: inout [Double], _ im: inout [Double], rows: vDSP_DFT_SetupD, columns: vDSP_DFT_SetupD) {
        precondition(re.count == nx * ny && im.count == nx * ny)
        let nx = self.nx, ny = self.ny
        re.withUnsafeMutableBufferPointer { r in
            im.withUnsafeMutableBufferPointer { m in
                outRe.withUnsafeMutableBufferPointer { orr in
                    outIm.withUnsafeMutableBufferPointer { oi in
                        // Along y: each row ix is contiguous.
                        for ix in 0..<nx {
                            let base = ix * ny
                            vDSP_DFT_ExecuteD(rows, r.baseAddress! + base, m.baseAddress! + base,
                                              orr.baseAddress!, oi.baseAddress!)
                            for iy in 0..<ny {
                                r[base + iy] = orr[iy]
                                m[base + iy] = oi[iy]
                            }
                        }
                        lineRe.withUnsafeMutableBufferPointer { lr in
                            lineIm.withUnsafeMutableBufferPointer { li in
                                // Along x: gather column iy.
                                for iy in 0..<ny {
                                    for ix in 0..<nx {
                                        lr[ix] = r[ix * ny + iy]
                                        li[ix] = m[ix * ny + iy]
                                    }
                                    vDSP_DFT_ExecuteD(columns, lr.baseAddress!, li.baseAddress!,
                                                      orr.baseAddress!, oi.baseAddress!)
                                    for ix in 0..<nx {
                                        r[ix * ny + iy] = orr[ix]
                                        m[ix * ny + iy] = oi[ix]
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
