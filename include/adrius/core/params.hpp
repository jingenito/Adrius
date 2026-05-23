// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <cstddef>
#include <cstdint>

namespace adrius {

struct GSOParams {
    // Squared norms below this value are treated as zero — signals linear
    // dependence. Lower values tolerate more floating-point noise; raise for
    // ill-conditioned inputs where false positives occur.
    double zero_threshold = 1e-14;
};

struct LLLParams {
    // Lovász condition parameter δ ∈ (0.25, 1.0). The classic LLL paper uses
    // 0.75; values closer to 1 produce a more reduced basis at more iterations.
    double delta = 0.75;

    // Size-reduction threshold η ∈ (0.5, 1.0).
    double eta = 0.51;

    // Maximum column swaps before throwing ConvergenceError; 0 = no limit.
    std::size_t max_iter = 0;

    GSOParams gso{};
};

struct ILLLParams {
    // Forwarded to each inner LLL call.
    LLLParams lll{};

    // Maximum outer ILLL iterations.
    std::size_t max_iterations = 30;

    // Stop once the denominator q exceeds this value; 0 = no limit.
    // The algorithm guarantees bounded Dirichlet coefficient, so useful
    // denominator ranges are typically 10^3–10^9.
    std::int64_t max_denominator = 1'000'000LL;

    // Bosma-Smeets ε ∈ (0, 1): controls the initial scale c₀ = ε^{n+1}
    // and the per-iteration decay 2^{−n(n+1)/4}.  Default ε = 0.5 (d = 2).
    double epsilon = 0.5;

    // Stop early when max_i |q·αᵢ − pᵢ| drops below this threshold.
    double quality_tol = 1e-12;
};

struct CFParams {
    // Maximum number of partial quotients to emit.
    std::size_t max_depth = 64;

    // Fractional parts below this threshold are treated as zero, triggering
    // early termination for inputs that are exact rationals.
    double zero_threshold = 1e-14;
};

struct RationalApproxParams {
    // The best-rational search will not exceed this denominator.
    std::int64_t max_denominator = 1'000'000LL;
};

} // namespace adrius
