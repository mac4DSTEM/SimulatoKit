// Diagnostic, step 2: is the round-trip loss a scalar gain (ψ_out = s ψ_in + small residual) and is s the same for different
// inputs (random, a probe-like smooth field, a single delta)? s fitted as Σ conj(in) out / Σ |in|².
import Foundation
import Metal
import MetalPerformanceShadersGraph
let device = MTLCreateSystemDefaultDevice()!
let queue = device.makeCommandQueue()!
let gd = MPSGraphDevice(mtlDevice: device)
func run(nx: Int, ny: Int, input a: [Float], label: String) {
  let g = MPSGraph(); let shape: [NSNumber] = [NSNumber(value: nx), NSNumber(value: ny)]
  let x = g.placeholder(shape: shape, dataType: .complexFloat32, name: nil)
  let f = MPSGraphFFTDescriptor(); f.inverse = false; f.scalingMode = .none
  let i = MPSGraphFFTDescriptor(); i.inverse = true; i.scalingMode = .size
  let y = g.fastFourierTransform(g.fastFourierTransform(x, axes: [0, 1], descriptor: f, name: nil), axes: [0, 1], descriptor: i, name: nil)
  let td = MPSGraphTensorData(device: gd, data: Data(bytes: a, count: a.count*4), shape: shape, dataType: .complexFloat32)
  let out = g.run(with: queue, feeds: [x: td], targetTensors: [y], targetOperations: nil)[y]!
  var b = a; out.mpsndarray().readBytes(&b, strideBytes: nil)
  var num_re = 0.0, num_im = 0.0, den = 0.0
  for q in stride(from: 0, to: a.count, by: 2) {
    let (ar, ai, br, bi) = (Double(a[q]), Double(a[q+1]), Double(b[q]), Double(b[q+1]))
    num_re += ar*br + ai*bi; num_im += ar*bi - ai*br; den += ar*ar + ai*ai
  }
  let sr = num_re/den, si = num_im/den
  var res = 0.0
  for q in stride(from: 0, to: a.count, by: 2) {
    let (ar, ai, br, bi) = (Double(a[q]), Double(a[q+1]), Double(b[q]), Double(b[q+1]))
    let pr = sr*ar - si*ai, pi = sr*ai + si*ar
    res += (br-pr)*(br-pr) + (bi-pi)*(bi-pi)
  }
  print(String(format: "%@ %dx%d: s-1 = %+.3e %+.3ei   residual energy / energy = %.3e   |s|^2-1 = %+.3e", label, nx, ny, sr-1, si, res/den, sr*sr+si*si-1))
}
for (nx, ny) in [(256, 256), (512, 512), (120, 160)] {
  var seed: UInt64 = 99
  var r = [Float](repeating: 0, count: 2*nx*ny)
  for k in r.indices { seed = seed &* 6364136223846793005 &+ 1442695040888963407; r[k] = Float(Double(seed >> 11) / Double(1 << 53)) - 0.5 }
  run(nx: nx, ny: ny, input: r, label: "random")
  var smooth = [Float](repeating: 0, count: 2*nx*ny)
  for ix in 0..<nx { for iy in 0..<ny { let dx = Double(ix - nx/2), dy = Double(iy - ny/3); let g = exp(-(dx*dx+dy*dy)/200); smooth[2*(ix*ny+iy)] = Float(g*cos(0.05*dx*dx)); smooth[2*(ix*ny+iy)+1] = Float(g*sin(0.05*dx*dx)) } }
  run(nx: nx, ny: ny, input: smooth, label: "smooth")
  var delta = [Float](repeating: 0, count: 2*nx*ny); delta[2*(7*ny+11)] = 1
  run(nx: nx, ny: ny, input: delta, label: "delta ")
  var ones = [Float](repeating: 0, count: 2*nx*ny); for k in stride(from: 0, to: ones.count, by: 2) { ones[k] = 1 }
  run(nx: nx, ny: ny, input: ones, label: "ones  ")
}
