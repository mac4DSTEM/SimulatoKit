// Gate D (registration amendment 2): is MPSGraph's float32 FFT round trip non-unitary by a systematic amount?
// Measures Σ|ψ|² after k round trips (forward .none, inverse .size) of a random complex field, for several sizes, and the
// same with forward .unitary / inverse .unitary, and with a unit-modulus multiply in between.
import Foundation
import Metal
import MetalPerformanceShadersGraph
let device = MTLCreateSystemDefaultDevice()!
let queue = device.makeCommandQueue()!
let gd = MPSGraphDevice(mtlDevice: device)
func experiment(nx: Int, ny: Int, mode: String) {
  let g = MPSGraph(); let shape: [NSNumber] = [NSNumber(value: nx), NSNumber(value: ny)]
  let x = g.placeholder(shape: shape, dataType: .complexFloat32, name: nil)
  let f = MPSGraphFFTDescriptor(); f.inverse = false; f.scalingMode = mode == "unitary" ? .unitary : .none
  let i = MPSGraphFFTDescriptor(); i.inverse = true; i.scalingMode = mode == "unitary" ? .unitary : .size
  let y = g.fastFourierTransform(g.fastFourierTransform(x, axes: [0, 1], descriptor: f, name: nil), axes: [0, 1], descriptor: i, name: nil)
  var a = [Float](repeating: 0, count: 2*nx*ny)
  var seed: UInt64 = 12345
  for k in a.indices { seed = seed &* 6364136223846793005 &+ 1442695040888963407; a[k] = Float(Double(seed >> 11) / Double(1 << 53)) - 0.5 }
  func energy(_ v: [Float]) -> Double { var s = 0.0; for q in stride(from: 0, to: v.count, by: 2) { s += Double(v[q])*Double(v[q]) + Double(v[q+1])*Double(v[q+1]) }; return s }
  let e0 = energy(a)
  var cur = MPSGraphTensorData(device: gd, data: Data(bytes: a, count: a.count*4), shape: shape, dataType: .complexFloat32)
  var line = "\(nx)x\(ny) \(mode):"
  for k in 1...100 {
    cur = g.run(with: queue, feeds: [x: cur], targetTensors: [y], targetOperations: nil)[y]!
    if [1, 10, 50, 100].contains(k) { var b = a; cur.mpsndarray().readBytes(&b, strideBytes: nil); line += String(format: "  k=%d dE/E=%+.3e", k, energy(b)/e0 - 1) }
  }
  print(line)
}
for (nx, ny) in [(128, 128), (256, 256), (512, 512), (120, 160)] { experiment(nx: nx, ny: ny, mode: "size"); experiment(nx: nx, ny: ny, mode: "unitary") }
