import Testing
import SimulatoDetector

/// The empty module builds and exposes its scope statement; replaced by real tests at its registration.
@Suite("SimulatoDetector placeholder")
struct SimulatoDetectorPlaceholderTests {
    @Test func plannedScopeIsStated() {
        #expect(!(SimulatoDetector.plannedScope.isEmpty))
    }
}
