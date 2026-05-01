#pragma once

#include <cstdint>
#include <functional>
#include <string_view>
#include <vector>

namespace bromesh {

/// One symbol in an L-system word, with optional parametric values.
struct Module {
    char symbol = 0;
    std::vector<float> params;
};

/// A single production rule. `predecessor` matches a symbol; `condition`
/// (if not null) further filters by parameters; `successor` produces the
/// replacement modules. Multiple rules with the same predecessor are
/// selected stochastically by `weight`.
struct ProductionRule {
    char predecessor = 0;
    std::function<bool(const std::vector<float>&)> condition;
    float weight = 1.0f;
    std::function<std::vector<Module>(const std::vector<float>&)> successor;
};

/// Parametric stochastic L-system rewriter. Produces module sequences
/// only — turtle interpretation and geometry building are caller concerns.
/// Bracket symbols `[` and `]` are passed through unchanged when no rule
/// matches them.
class LSystem {
public:
    void addRule(ProductionRule rule);
    void setAxiom(std::vector<Module> axiom);

    /// Run `iterations` rewrite passes. Deterministic given `seed`.
    std::vector<Module> derive(int iterations, uint64_t seed = 0) const;

    const std::vector<Module>& axiom() const { return axiom_; }

private:
    std::vector<Module> axiom_;
    std::vector<ProductionRule> rules_;
};

/// Parse a compact textual form like `F(1.0)[+(25)F]F` into a module list.
/// A symbol is any non-whitespace char that isn't `(`, `)`, or `,`. Each
/// symbol may be followed by a parenthesized comma-separated list of
/// floats. Whitespace is ignored. Returns an empty vector on parse error.
std::vector<Module> parseModules(std::string_view s);

} // namespace bromesh
