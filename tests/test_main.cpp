#include "test_framework.h"

// Entry point + counters. Tests register themselves via static constructors
// in test_mesh.cpp / test_animation.cpp / test_rigging.cpp, but run from
// main() — not during static init — so cross-TU init-order is not a factor.

int tests_run = 0;
int tests_passed = 0;

std::vector<TestEntry>& testRegistry() {
    static std::vector<TestEntry> v;
    return v;
}

int main() {
    for (const auto& t : testRegistry()) t.fn();
    std::printf("bromesh tests: %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
