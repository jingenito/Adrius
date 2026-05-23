// Copyright (c) 2025 InGenifold Research LLC. MIT License.

#include <adrius/linalg/lll.hpp>

#include <gtest/gtest.h>
#include <cmath>

using B = adrius::EigenBackend;

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

TEST(LLL, AlreadyReducedBasisZeroSwaps) {
    // Identity is LLL-reduced for any δ; expect zero swaps.
    auto basis = make_basis({{1,0,0},{0,1,0},{0,0,1}});
    auto result = adrius::lll_reduce<B>(basis);
    EXPECT_EQ(result.swap_count, 0u);
}

TEST(LLL, TransformIsUnimodular) {
    // |det(U)| must be 1 for any unimodular matrix.
    auto basis = make_basis({{1,1,1},{0,1,1},{0,0,1}});
    auto result = adrius::lll_reduce<B>(basis);

    // Convert int_matrix to double for determinant via Eigen
    const std::size_t n = B::cols(result.reduced_basis);
    adrius::EigenBackend::matrix_type Ud =
        adrius::EigenBackend::matrix_type::Zero(
            static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(n));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            Ud(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
                static_cast<double>(B::int_get(result.transform, i, j));

    EXPECT_NEAR(std::abs(Ud.determinant()), 1.0, 1e-9);
}

TEST(LLL, ReducedBasisEqualsOriginalTimesU) {
    // reduced_basis = original_basis · U must hold exactly.
    auto original = make_basis({{2,1},{1,2}});
    auto result   = adrius::lll_reduce<B>(original);

    const std::size_t m = B::rows(original);
    const std::size_t n = B::cols(original);

    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double reconstructed = 0.0;
            for (std::size_t k = 0; k < n; ++k)
                reconstructed += B::get(original, i, k)
                                 * static_cast<double>(B::int_get(result.transform, k, j));
            EXPECT_NEAR(reconstructed, B::get(result.reduced_basis, i, j), 1e-9)
                << "Mismatch at (" << i << ", " << j << ")";
        }
    }
}

TEST(LLL, ClassicLovaszExample) {
    // From Lenstra, Lenstra & Lovász (1982). The input basis:
    //   b0 = [1, 1, 1], b1 = [-1, 0, 2], b2 = [3, 5, 6]
    // LLL-reduces to shorter vectors. Just verify swap_count > 0 and
    // the Lovász condition holds on the output.
    auto basis = make_basis({{1,-1,3},{1,0,5},{1,2,6}});
    adrius::LLLParams params{.delta = 0.75};
    auto result = adrius::lll_reduce<B>(basis, params);

    // Every consecutive pair must satisfy the Lovász condition.
    auto gso = adrius::gram_schmidt<B>(result.reduced_basis);
    const std::size_t n = B::cols(result.reduced_basis);
    for (std::size_t k = 1; k < n; ++k) {
        const double mu_k = B::get(gso.mu, k, k - 1);
        const double rhs  = (params.delta - mu_k * mu_k) * gso.norms_sq[k - 1];
        EXPECT_GE(gso.norms_sq[k], rhs - 1e-9)
            << "Lovász condition violated at k=" << k;
    }
}
