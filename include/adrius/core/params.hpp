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

    // Maximum outer ILLL iterations; 0 = no limit.
    std::size_t max_iterations = 0;

    // Terminate early when the best approximation quality (||q·α - p||_∞ · q)
    // stops improving by more than this relative fraction between iterations.
    double convergence_tol = 1e-10;
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
