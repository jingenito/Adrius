// Copyright (c) 2025 InGenifold Research LLC. MIT License.

#include <adrius/linalg/gram_schmidt.hpp>

#include <gtest/gtest.h>
#include <cmath>

using B = adrius::EigenBackend;

// Build a 3×3 basis from column vectors given as initializer lists.
static B::matrix_type make_basis(std::initializer_list<std::initializer_list<double>> cols) {
    const std::size_t m = cols.begin()->size();
    const std::size_t n = cols.size();
    auto M = B::make_zero_matrix(m, n);
    std::size_t j = 0;
    for (auto& col : cols) {
        std::size_t i = 0;
        for (double v : col) { B::set(M, i++, j, v); }
        ++j;
    }
    return M;
}

TEST(GramSchmidt, OrthogonalInputIsUnchanged) {
    // Identity columns are already orthogonal; GSO should be identity.
    auto basis = make_basis({{1,0,0},{0,1,0},{0,0,1}});
    auto result = adrius::gram_schmidt<B>(basis);

    ASSERT_EQ(result.norms_sq.size(), 3u);
    for (double ns : result.norms_sq)
        EXPECT_DOUBLE_EQ(ns, 1.0);

    // mu should be all zeros (no projection coefficients)
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < i; ++j)
            EXPECT_DOUBLE_EQ(B::get(result.mu, i, j), 0.0);
}

TEST(GramSchmidt, NonOrthogonalBasis) {
    // b0 = [1,1,0], b1 = [1,0,0], b2 = [0,0,1]
    // b*0 = [1,1,0],  ||b*0||^2 = 2
    // mu(1,0) = <[1,0,0],[1,1,0]>/2 = 0.5
    // b*1 = [1,0,0] - 0.5*[1,1,0] = [0.5,-0.5,0], ||b*1||^2 = 0.5
    // b*2 = [0,0,1] (already orthogonal to b*0 and b*1)
    auto basis = make_basis({{1,1,0},{1,0,0},{0,0,1}});
    auto result = adrius::gram_schmidt<B>(basis);

    EXPECT_DOUBLE_EQ(result.norms_sq[0], 2.0);
    EXPECT_DOUBLE_EQ(result.norms_sq[1], 0.5);
    EXPECT_DOUBLE_EQ(result.norms_sq[2], 1.0);

    EXPECT_DOUBLE_EQ(B::get(result.mu, 1, 0), 0.5);
    EXPECT_DOUBLE_EQ(B::get(result.mu, 2, 0), 0.0);
    EXPECT_DOUBLE_EQ(B::get(result.mu, 2, 1), 0.0);
}

TEST(GramSchmidt, OrthogonalityOfOutput) {
    // Any two distinct output columns must be orthogonal.
    auto basis = make_basis({{3,1,0},{2,2,0},{0,0,5}});
    auto result = adrius::gram_schmidt<B>(basis);

    const std::size_t n = 3;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j)
            EXPECT_NEAR(B::inner_product(B::col(result.Q, i), B::col(result.Q, j)),
                        0.0, 1e-12)
                << "Columns " << i << " and " << j << " are not orthogonal";
}

TEST(GramSchmidt, LinearDependenceThrows) {
    // Columns 0 and 1 are the same — must throw DomainError.
    auto basis = make_basis({{1,0,0},{1,0,0},{0,1,0}});
    // Cast to void inside the macro to silence [[nodiscard]] on MSVC (/W4).
    // The exception is thrown before the cast executes, so it is still caught.
    EXPECT_THROW((void)adrius::gram_schmidt<B>(basis), adrius::DomainError);
}
