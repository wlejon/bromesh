#include "bromesh/procedural/lsystem.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <unordered_map>

namespace bromesh {

void LSystem::addRule(ProductionRule rule) {
    rules_.push_back(std::move(rule));
}

void LSystem::setAxiom(std::vector<Module> axiom) {
    axiom_ = std::move(axiom);
}

std::vector<Module> LSystem::derive(int iterations, uint64_t seed) const {
    std::vector<Module> out;
    deriveWithin(iterations, seed, SIZE_MAX, out);
    return out;
}

bool LSystem::deriveWithin(int iterations, uint64_t seed, size_t maxModules,
                           std::vector<Module>& out, int* passes) const {
    if (passes) *passes = 0;
    out.clear();
    if (axiom_.size() > maxModules) return false;
    std::vector<Module> current = axiom_;
    if (iterations <= 0 || rules_.empty()) {
        out = std::move(current);
        return true;
    }

    std::mt19937_64 rng(seed);
    std::vector<Module> next;

    for (int it = 0; it < iterations; ++it) {
        next.clear();
        next.reserve(std::min(current.size() * 2, maxModules));
        for (const Module& m : current) {
            // Collect matching rules.
            float total = 0.0f;
            std::vector<const ProductionRule*> matches;
            matches.reserve(4);
            for (const ProductionRule& r : rules_) {
                if (r.predecessor != m.symbol) continue;
                if (r.condition && !r.condition(m.params)) continue;
                if (!r.successor) continue;
                matches.push_back(&r);
                total += (r.weight > 0.0f) ? r.weight : 0.0f;
            }
            if (matches.empty() || total <= 0.0f) {
                if (next.size() >= maxModules) return false;
                next.push_back(m);
                continue;
            }
            // Weighted pick.
            std::uniform_real_distribution<float> dist(0.0f, total);
            float pick = dist(rng);
            const ProductionRule* chosen = matches.front();
            float acc = 0.0f;
            for (const ProductionRule* r : matches) {
                acc += (r->weight > 0.0f) ? r->weight : 0.0f;
                if (pick <= acc) {
                    chosen = r;
                    break;
                }
            }
            std::vector<Module> rep = chosen->successor(m.params);
            if (rep.size() > maxModules - next.size()) return false;
            for (Module& rm : rep) next.push_back(std::move(rm));
        }
        current.swap(next);
        if (passes) *passes = it + 1;
    }
    out = std::move(current);
    return true;
}

std::vector<Module> parseModules(std::string_view s) {
    std::vector<Module> out;
    size_t i = 0;
    auto skipWs = [&]() {
        while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    };
    while (i < s.size()) {
        skipWs();
        if (i >= s.size()) break;
        char c = s[i++];
        if (c == '(' || c == ')' || c == ',') {
            // Stray punctuation outside a parameter list — parse error.
            return {};
        }
        Module m;
        m.symbol = c;
        skipWs();
        if (i < s.size() && s[i] == '(') {
            ++i;
            while (true) {
                skipWs();
                if (i >= s.size()) return {};
                if (s[i] == ')') { ++i; break; }
                // Parse a float.
                size_t start = i;
                while (i < s.size() && s[i] != ',' && s[i] != ')'
                       && !std::isspace(static_cast<unsigned char>(s[i]))) {
                    ++i;
                }
                if (start == i) return {};
                std::string token(s.substr(start, i - start));
                char* end = nullptr;
                float v = std::strtof(token.c_str(), &end);
                if (end == token.c_str()) return {};
                m.params.push_back(v);
                skipWs();
                if (i < s.size() && s[i] == ',') { ++i; continue; }
            }
        }
        out.push_back(std::move(m));
    }
    return out;
}

std::string serializeModules(const std::vector<Module>& mods) {
    std::string s;
    s.reserve(mods.size() * 2);
    char buf[64];
    for (const auto& m : mods) {
        s.push_back(m.symbol);
        if (!m.params.empty()) {
            s.push_back('(');
            for (size_t i = 0; i < m.params.size(); i++) {
                if (i) s.push_back(',');
                std::snprintf(buf, sizeof(buf), "%g", static_cast<double>(m.params[i]));
                s += buf;
            }
            s.push_back(')');
        }
    }
    return s;
}

} // namespace bromesh

