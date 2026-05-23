// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/core/concepts.hpp>
#include <adrius/core/error.hpp>
#include <adrius/core/params.hpp>
#include <adrius/core/traits.hpp>
#include <adrius/backend/default_backend.hpp>

#include <string>
#include <vector>

namespace adrius {

// Result of Gram-Schmidt orthogonalization.
//
//   Q:        columns are the orthogonalized vectors b*_i (not normalized)
//   mu:       GSO coefficient matrix (lower-triangular, n×n):
//               mu(i, j) = <b_i, b*_j> / ||b*_j||²   for j < i
//   norms_sq: squared norms ||b*_i||² in index order
//
// LLL stores mu and norms_sq directly and updates them incrementally
// using the Lovász swap formulas — never performing a full recompute.
template <Backend B>
struct [[nodiscard]] GSOResult {
    matrix_of<B>              Q;
    matrix_of<B>              mu;
    std::vector<scalar_of<B>> norms_sq;

    GSOResult() = default;
    GSOResult(GSOResult&&) noexcept = default;
    GSOResult& operator=(GSOResult&&) noexcept = default;
    GSOResult(const GSOResult&) = default;
    GSOResult& operator=(const GSOResult&) = default;
};

// Compute classical Gram-Schmidt orthogonalization on the columns of `basis`.
//
// Throws DomainError when a column projects to a near-zero vector, which
// indicates the basis is linearly dependent.
template <Backend B = DefaultBackend>
[[nodiscard]] GSOResult<B>
gram_schmidt(const matrix_of<B>& basis, GSOParams params = {})
{
    const std::size_t n = B::cols(basis);
    const std::size_t m = B::rows(basis);
    const scalar_of<B> threshold = static_cast<scalar_of<B>>(params.zero_threshold);

    GSOResult<B> result;
    result.Q        = B::make_zero_matrix(m, n);
    result.mu       = B::make_zero_matrix(n, n);
    result.norms_sq.resize(n, scalar_of<B>{0});

    for (std::size_t i = 0; i < n; ++i) {
        // b*_i = b_i - Σ_{j<i} mu(i,j)·b*_j
        vector_of<B> qi = B::col(basis, i);

        for (std::size_t j = 0; j < i; ++j) {
            const scalar_of<B> mu_ij =
                B::inner_product(B::col(basis, i), B::col(result.Q, j))
                / result.norms_sq[j];
            B::set(result.mu, i, j, mu_ij);
            qi = B::subtract(qi, B::scale(B::col(result.Q, j), mu_ij));
        }

        const scalar_of<B> ns = B::norm_sq(qi);
        if (ns < threshold)
            throw DomainError(
                "gram_schmidt: basis is linearly dependent at column "
                + std::to_string(i));

        result.norms_sq[i] = ns;
        B::set_col(result.Q, i, qi);
    }

    return result;
}

} // namespace adrius
