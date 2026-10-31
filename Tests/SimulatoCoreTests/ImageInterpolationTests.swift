//
//  ImageInterpolationTests.swift
//  Role: Fourier-space image interpolation (abTEM's Images.interpolate, method "fft") and integrate_disc(interpolate=),
//        against numbers abTEM 1.0.10 printed in float64 (tools/pins/image-interpolate/pins.py) for a planted 5 x 6 image
//        at sampling (0.3, 0.25) A, values sin(0.7 i) + cos(0.4 j)^2 + 0.05 i j; and the refusals.
//

import Foundation
import Testing
import SimulatoCore

@Suite("Image interpolation")
struct ImageInterpolationTests {
    static let image: RealField = {
        var v = [Double]()
        for i in 0..<5 { for j in 0..<6 { v.append(sin(0.7 * Double(i)) + pow(cos(0.4 * Double(j)), 2) + 0.05 * Double(i * j)) } }
        return RealField(nx: 5, ny: 6, values: v, samplingX: 0.3, samplingY: 0.25)
    }()

    /// assert_allclose defaults, rtol 1e-7 and atol 0 (abTEM states no tolerance for interpolate or integrate_disc); a
    /// value that is zero in exact arithmetic is held at 1e-12 absolute.
    static func close(_ a: [Double], _ b: [Double], rtol: Double = 1e-7) -> Bool {
        a.count == b.count && zip(a, b).allSatisfy { abs($0 - $1) <= rtol * abs($1) + 1e-12 }
    }

    static func samplingIs(_ f: RealField, _ x: Double, _ y: Double) -> Bool {
        abs(f.samplingX - x) <= 1e-12 * x && abs(f.samplingY - y) <= 1e-12 * y
    }

    @Test func upsamplingBySamplingMatchesAbTEM() throws {
        let f = try Self.image.interpolated(samplingX: 0.2, samplingY: 0.2)
        #expect(f.nx == 8 && f.ny == 8 && Self.samplingIs(f, 0.1875, 0.1875))
        #expect(Self.close(f.values, [
            1.0, 0.9601480812272941, 0.6488274077236056, 0.39300261176602336, 0.13130314222937708, 0.02960662516674309,
            -0.020738858337560846, 0.37663437345536155, 1.3578464522963176, 1.3179037027850447, 1.0062075221487785,
            0.7502671381805395, 0.4883582744926171, 0.3865296207025732, 0.33572751841938686, 0.7334393056889427,
            1.8003221607478836, 1.7884250055000486, 1.5926732820974674, 1.3724227372710842, 1.1751680149065848,
            1.1141389087442517, 1.2043259166194342, 1.4975027228559856, 1.980309754460406, 1.9759049882857542,
            1.8111278620449929, 1.6004118686457733, 1.4204295198825416, 1.3703000253299793, 1.4981523170449673,
            1.7634025711714012, 1.9544043049917348, 1.9489875331015802, 1.7800266319584528, 1.5680227955968893,
            1.385707447221112, 1.3341057294380017, 1.4568705274110618, 1.7258928534458902, 1.8257019752052555,
            1.8436225234714048, 1.771141407830237, 1.5888358261947761, 1.4603205285978382, 1.44266895985691,
            1.6827537180137984, 1.864790320124038, 1.507565481485125, 1.5430642398352024, 1.5432539333120896,
            1.3833177650113644, 1.295325911785178, 1.303246381660814, 1.6316992512091162, 1.7482161912841276,
            1.0984337652627647, 1.1031936451219218, 0.9763044630676065, 0.7772510766242722, 0.618396193169255,
            0.5815990391203479, 0.7555232834238985, 0.9866139440569834]))
    }

    @Test func croppingAndMixedResamplingMatchAbTEM() throws {
        let c = try Self.image.interpolated(rows: 3, columns: 4)
        #expect(Self.samplingIs(c, 0.5, 0.375))
        #expect(Self.close(c.values, [1.0198993161247643, 0.8634975257608967, 0.3623418465144165, 0.36713634045661836,
                                      1.8778708201312442, 1.7255956923912834, 1.2511151796771314, 1.282584840151661,
                                      1.7435342696058869, 1.6231981954827803, 1.355174970839304, 1.5931019193845097]))
        let m = try Self.image.interpolated(samplingX: 0.17, samplingY: 0.4)
        #expect(m.nx == 9 && m.ny == 4 && Self.samplingIs(m, 0.16666666666666666, 0.375))
        #expect(Self.close(m.values, [
            0.8815284818144458, 0.7117292318290219, 0.12397101220409798, 0.042162965767856046, 1.1788714665239841,
            1.0071376595233525, 0.40687426677184824, 0.31256104720902594, 1.6653656652749362, 1.512975660804376,
            1.037752573230411, 1.0684796588451277, 1.9384347446984385, 1.7966282406645506, 1.3898179753570876,
            1.4889578832383075, 1.9566389495907732, 1.8124523223885252, 1.3902566990642402, 1.4740112489286374,
            1.9163835149498045, 1.779214369564729, 1.4023804617890283, 1.5314967272020097, 1.8213411793490102,
            1.7039339411413879, 1.4548430094696656, 1.711702250986625, 1.5057939897471375, 1.3927014317230826,
            1.1715010311947631, 1.4562508038551254, 1.0595552256371534, 0.9201013832658551, 0.5284989620114121,
            0.6428467139456514]))
    }

    @Test func intensityNormalizationKeepsTheSumAndMatchesAbTEM() throws {
        let f = try Self.image.interpolated(rows: 9, columns: 7, normalization: .intensity)
        #expect(Self.samplingIs(f, 0.16666666666666666, 0.21428571428571427))
        #expect(Self.close(f.values, [
            0.47619047619047616, 0.4349265773060213, 0.2758936847900825, 0.1270881728466438, 0.029723204763476653,
            -0.012213978011907439, 0.1345516051280584, 0.6212204046573486, 0.5777820583619322, 0.41448203892473934,
            0.2632183404911704, 0.1637897420018617, 0.11648576161200819, 0.26422535397444846, 0.8185071403226568,
            0.7968112813248561, 0.6761786321882627, 0.5494945525890664, 0.4707003726036572, 0.47705946466763216,
            0.6150598547662446, 0.9297313870235687, 0.9199314116531432, 0.8226431970715469, 0.7094072699885073,
            0.6419027202041407, 0.6776222831406525, 0.8102941008837611, 0.942629932492581, 0.9301546921283996,
            0.8276165481997899, 0.7113562664085074, 0.6413127919494073, 0.670429479270632, 0.8042996395274963,
            0.9109894404336063, 0.9064018689511728, 0.8193424549571537, 0.7119990874062382, 0.6494413002329724,
            0.6980257025814224, 0.8283627066920775, 0.8306109750805063, 0.8482358416703951, 0.8047660275017299,
            0.7225335524502049, 0.6810561843155885, 0.7844635541256281, 0.9048508489781881, 0.6726825011446215,
            0.695157080139228, 0.6612043222388302, 0.5844543883856783, 0.5475795753319368, 0.662956618371733,
            0.7811715614880629, 0.5070362575382878, 0.4998806803535306, 0.407781822217943, 0.29753535529005154,
            0.23254043644657624, 0.27478669200531863, 0.40627399353168553]))
        #expect(abs(f.sum - Self.image.sum) <= 1e-12 * Self.image.sum)
        // The same size twice is the identity.
        #expect(Self.close(try Self.image.interpolated(rows: 5, columns: 6).values, Self.image.values))
    }

    @Test func discOfAnInterpolatedImageMatchesAbTEM() throws {
        let p = Self.image
        #expect(Self.close([try p.integrateDisc(x: 0.7, y: 0.8, radius: 0.45, interpolate: 0.1)], [84.58187215553677]))
        #expect(Self.close([try p.integrateDisc(x: 0.7, y: 0.8, radius: 0.45, mean: true, interpolate: 0.1)],
                           [1.4080820584300002]))
        #expect(Self.close([try p.integrateDisc(x: 0.1, y: 1.4, radius: 0.5, interpolate: 0.1)], [75.40052718074456]))
        #expect(Self.close([try p.integrateDisc(x: 0.7, y: 0.8, radius: 0.45)], [10.906912550258145]))
    }

    @Test func refusals() {
        let p = Self.image
        #expect(throws: SimulatoError.self) { try p.interpolated(rows: 1, columns: 4) }
        #expect(throws: SimulatoError.self) { try p.interpolated(rows: 4, columns: 0) }
        #expect(throws: SimulatoError.self) { try p.interpolated(samplingX: 0, samplingY: 0.1) }
        #expect(throws: SimulatoError.self) { try p.interpolated(samplingX: 0.1, samplingY: .nan) }
        #expect(throws: SimulatoError.self) { try p.interpolated(samplingX: 1e-300, samplingY: 0.1) }
        #expect(throws: SimulatoError.self) { try p.interpolated(rows: 1 << 40, columns: 1 << 40) }
        #expect(throws: SimulatoError.self) { try p.integrateDisc(x: 0.7, y: 0.8, radius: 0.45, interpolate: -1) }
        // A one-pixel axis is refused (abTEM keeps its zero frequency; the kit's FFT starts at 2 points).
        let line = RealField(nx: 1, ny: 4, values: [1, 2, 3, 4], samplingX: 0.3, samplingY: 0.25)
        #expect(throws: SimulatoError.self) { try line.interpolated(rows: 2, columns: 8) }
    }
}
