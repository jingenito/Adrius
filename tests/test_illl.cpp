// Copyright (c) 2025 InGenifold Research LLC. MIT License.

// ILLL tests are stubs until the Bosma-Smeets scale update rule (Theorem 23,
// eq. 24 / Lemma 25, eq. 28) is implemented. Each test documents what the
// expected invariant is so the implementation target is unambiguous.

#include <adrius/approx/illl.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <numbers>

using B = adrius::EigenBackend;

TEST(ILLL, PreprocessProducesCorrectLatticeShape) {
    // For n inputs, the lattice must be (n+1)×(n+1).
    const std::vector<double> alpha = {std::numbers::sqrt2, std::numbers::pi};
    auto prepared = adrius::preprocess_illl<B>(
        std::span<const double>{alpha}, 1e6);

    EXPECT_EQ(B::rows(prepared.lattice), 3u);
    EXPECT_EQ(B::cols(prepared.lattice), 3u);
}

TEST(ILLL, PreprocessFirstColumnIsScaledAlpha) {
    // First column of the lattice must be (N, ⌊N·α_1⌋, ⌊N·α_2⌋)ᵀ.
    const double N = 1000.0;
    const std::vector<double> alpha = {0.5, 0.25};
    auto prepared = adrius::preprocess_illl<B>(
        std::span<const double>{alpha}, N);

    EXPECT_DOUBLE_EQ(B::get(prepared.lattice, 0, 0), N);
    EXPECT_NEAR(B::get(prepared.lattice, 1, 0), std::round(N * alpha[0]), 1.0);
    EXPECT_NEAR(B::get(prepared.lattice, 2, 0), std::round(N * alpha[1]), 1.0);
}

TEST(ILLL, PreprocessIdentityBlockIsCorrect) {
    // Columns 1..n of the lattice must form an identity block.
    const std::vector<double> alpha = {0.3, 0.7};
    auto prepared = adrius::preprocess_illl<B>(
        std::span<const double>{alpha}, 1e4);

    for (std::size_t j = 1; j <= 2; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            EXPECT_DOUBLE_EQ(B::get(prepared.lattice, i, j),
                             (i == j ? 1.0 : 0.0))
                << "Identity block mismatch at (" << i << "," << j << ")";
}

// TODO: add test once scale update rule is implemented
// TEST(ILLL, KnownSimultaneousApproximation) {
//     // For alpha = {sqrt(2)}, the best approximation with denominator <= 1000
//     // should be 1393/985 (a convergent of sqrt(2)).
//     // Expected: result.relations gives a relation with |q*sqrt(2) - p| < 1e-3.
// }
