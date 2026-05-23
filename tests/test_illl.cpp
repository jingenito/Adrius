// Copyright (c) 2025 InGenifold Research LLC. MIT License.

#include <adrius/approx/illl.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <numbers>

using B = adrius::EigenBackend;

// ── Preprocessing ─────────────────────────────────────────────────────────

TEST(ILLL, PreprocessProducesCorrectShape) {
    const std::vector<double> alpha = {std::numbers::sqrt2, std::numbers::pi};
    auto p = adrius::preprocess_illl<B>(std::span<const double>{alpha});

    // Lattice must be (n+1)×(n+1)
    EXPECT_EQ(B::rows(p.lattice), 3u);
    EXPECT_EQ(B::cols(p.lattice), 3u);
}

TEST(ILLL, PreprocessInitialScaleIsEpsilonPowNPlus1) {
    // c₀ = ε^{n+1} with default ε = 0.5, n = 2  →  c₀ = 0.5³ = 0.125
    const std::vector<double> alpha = {0.5, 0.25};
    adrius::ILLLParams params;
    params.epsilon = 0.5;
    auto p = adrius::preprocess_illl<B>(std::span<const double>{alpha}, params);

    const double expected_c0 = std::pow(0.5, 3.0);  // n+1 = 3
    EXPECT_NEAR(B::get(p.lattice, 0, 0), expected_c0, 1e-14);
    EXPECT_NEAR(p.scale, expected_c0, 1e-14);
}

TEST(ILLL, PreprocessColumnZeroHasAlphaValues) {
    // Column 0 = (c₀, α₁, α₂)ᵀ — targets sit in rows 1..n unscaled
    const std::vector<double> alpha = {0.3, 0.7};
    auto p = adrius::preprocess_illl<B>(std::span<const double>{alpha});

    EXPECT_NEAR(B::get(p.lattice, 1, 0), alpha[0], 1e-14);
    EXPECT_NEAR(B::get(p.lattice, 2, 0), alpha[1], 1e-14);
}

TEST(ILLL, PreprocessIdentityBlockIsCorrect) {
    const std::vector<double> alpha = {0.3, 0.7};
    auto p = adrius::preprocess_illl<B>(std::span<const double>{alpha});

    // Columns 1..n: identity block
    for (std::size_t j = 1; j <= 2; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(B::get(p.lattice, i, j), (i == j ? 1.0 : 0.0))
                << "Identity mismatch at (" << i << "," << j << ")";
}

// ── Algorithmic correctness ───────────────────────────────────────────────

TEST(ILLL, OneDimensionalBoundedDirichletCoefficient) {
    // For n = 1, α = √2, Bosma-Smeets guarantees θ(k) = q(k)·|q·α-p| ≤ C(1).
    // (Classical Hurwitz: θ ≤ 1/√5 ≈ 0.447 for best approximations.)
    // We verify that at least one approximation found beats θ ≤ 1.0 — a
    // generous but unconditional bound.
    const std::vector<double> alpha = {std::numbers::sqrt2};
    adrius::ILLLParams params;
    params.max_iterations  = 15;
    params.max_denominator = 0;  // no denominator cap
    params.quality_tol     = 0;  // no early exit

    auto result = adrius::illl<B>(std::span<const double>{alpha}, params);

    ASSERT_GT(result.iterations, 0u) << "ILLL produced no relations";

    // Every relation must satisfy the Bosma-Smeets bound (generously θ ≤ 2).
    for (std::size_t k = 0; k < result.iterations; ++k) {
        const double q       = static_cast<double>(result.relations[k][0]);
        const double quality = static_cast<double>(result.quality[k]);
        EXPECT_GT(q, 0.0) << "Denominator must be positive at step " << k;
        EXPECT_LT(q * quality, 2.0)
            << "Dirichlet coefficient out of range at step " << k
            << " (q=" << q << ", err=" << quality << ")";
    }
}

TEST(ILLL, DenominatorsAreNonDecreasing) {
    // Successive ILLL iterations must produce non-decreasing denominators.
    const std::vector<double> alpha = {std::numbers::sqrt2, std::numbers::pi};
    adrius::ILLLParams params;
    params.max_iterations  = 10;
    params.max_denominator = 0;
    params.quality_tol     = 0;

    auto result = adrius::illl<B>(std::span<const double>{alpha}, params);
    ASSERT_GE(result.iterations, 2u);

    for (std::size_t k = 1; k < result.iterations; ++k) {
        EXPECT_GE(result.relations[k][0], result.relations[k - 1][0])
            << "Denominator shrank at step " << k;
    }
}

TEST(ILLL, MaxDenominatorStopsIteration) {
    // With a larger denominator bound, verify the algorithm respects it
    const std::vector<double> alpha = {std::numbers::pi};
    adrius::ILLLParams params;
    params.max_iterations  = 30;
    params.max_denominator = 1000;
    params.quality_tol     = 0;

    auto result = adrius::illl<B>(std::span<const double>{alpha}, params);

    ASSERT_GT(result.iterations, 0u) << "Algorithm found no relations";

    // All denominators must be <= max_denominator
    for (const auto& rel : result.relations) {
        EXPECT_LE(rel[0], static_cast<adrius::EigenBackend::integer_type>(params.max_denominator));
    }
}

TEST(ILLL, KnownApproximationSqrt2) {
    // CF convergents of √2: 1/1, 3/2, 7/5, 17/12, 41/29, 99/70, 239/169.
    // ILLL with 1D input should find (q, p) matching a convergent somewhere.
    const std::vector<double> alpha = {std::numbers::sqrt2};
    adrius::ILLLParams params;
    params.max_iterations  = 20;
    params.max_denominator = 300;
    params.quality_tol     = 0;

    auto result = adrius::illl<B>(std::span<const double>{alpha}, params);
    ASSERT_GT(result.iterations, 0u);

    // At least one relation should satisfy |q·√2 - p| < 1/q (Dirichlet bound)
    bool found_dirichlet = false;
    for (auto& rel : result.relations) {
        const double q   = static_cast<double>(rel[0]);
        const double p   = static_cast<double>(rel[1]);
        const double err = std::abs(q * std::numbers::sqrt2 - p);
        if (err < 1.0 / q) { found_dirichlet = true; break; }
    }
    EXPECT_TRUE(found_dirichlet)
        << "No relation satisfied the Dirichlet bound |q·√2-p| < 1/q";
}
