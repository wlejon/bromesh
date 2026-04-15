#include "test_framework.h"

// Entry point + counters. Tests register themselves via static constructors
// scattered across test_mesh.cpp / test_animation.cpp / test_rigging.cpp and
// run before main() executes.

int tests_run = 0;
int tests_passed = 0;

int main() {
    std::printf("bromesh tests: %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
