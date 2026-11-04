#include "test_framework.h"

#include <cstring>

// Entry point + counters. Tests register themselves via static constructors
// in test_mesh.cpp / test_animation.cpp / test_rigging.cpp, but run from
// main() — not during static init — so cross-TU init-order is not a factor.

int tests_run = 0;
int tests_passed = 0;

std::vector<TestEntry>& testRegistry() {
    static std::vector<TestEntry> v;
    return v;
}

int main(int argc, char** argv) {
    const char* filter = nullptr;
    bool listOnly = false;

    for (int i = 1; i < argc; ++i) {
        const char* arg = argv[i];
        if (std::strcmp(arg, "--list") == 0 || std::strcmp(arg, "-l") == 0) {
            listOnly = true;
        } else if (std::strncmp(arg, "--filter=", 9) == 0) {
            filter = arg + 9;
        } else if ((std::strcmp(arg, "--filter") == 0 || std::strcmp(arg, "-f") == 0) && i + 1 < argc) {
            filter = argv[++i];
        } else if (arg[0] != '-') {
            filter = arg;
        }
    }

    if (listOnly) {
        for (const auto& t : testRegistry()) {
            if (filter && std::strstr(t.name, filter) == nullptr) continue;
            std::printf("%s\n", t.name);
        }
        return 0;
    }

    for (const auto& t : testRegistry()) {
        if (filter && std::strstr(t.name, filter) == nullptr) continue;
        t.fn();
    }
    std::printf("bromesh tests: %d/%d passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
