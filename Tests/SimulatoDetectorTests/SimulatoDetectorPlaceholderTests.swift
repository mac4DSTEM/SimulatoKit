import XCTest
import SimulatoDetector

/// The empty module builds and exposes its scope statement; replaced by real tests at its registration.
final class SimulatoDetectorPlaceholderTests: XCTestCase {
    func testPlannedScopeIsStated() {
        XCTAssertFalse(SimulatoDetector.plannedScope.isEmpty)
    }
}
