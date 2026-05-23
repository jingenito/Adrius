// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/core/concepts.hpp>
#include <adrius/core/error.hpp>
#include <adrius/core/params.hpp>
#include <adrius/core/traits.hpp>
#include <adrius/backend/default_backend.hpp>
#include <adrius/linalg/gram_schmidt.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace adrius {

// Output of lll_reduce.
//
//   reduced_basis: the LLL-reduced basis (columns)
//   transform:     unimodular integer matrix U such that reduced_basis = original_basis · U
//   swap_count:    number of Lovász swaps performed (useful for profiling)
template <Backend B>
struct [[nodiscard]] LLLResult {
    matrix_of<B>     reduced_basis;
    int_matrix_of<B> transform;
    std::size_t      swap_count{0};
    // Compiler-generated move/copy are correct. Explicit = default declarations
    // prevent aggregate initialization on MSVC, so we omit them here.
};

// Pre-computed state for lll_reduce: holds the working basis copy and its
// initial Gram-Schmidt decomposition. Move-only — copying a large basis is
// almost never what the caller wants.
template <Backend B>
struct PreparedLLL {
    matrix_of<B>     basis;
    int_matrix_of<B> transform;  // starts as identity
    GSOResult<B>     gso;

    PreparedLLL() = default;
    PreparedLLL(PreparedLLL&&) noexcept = default;
    PreparedLLL& operator=(PreparedLLL&&) noexcept = default;
    PreparedLLL(const PreparedLLL&) = delete;
    PreparedLLL& operator=(const PreparedLLL&) = delete;
};

// Stage 1: validate the basis and compute the initial GSO.
// Throws DomainError on a linearly dependent basis.
template <Backend B = DefaultBackend>
[[nodiscard]] PreparedLLL<B>
preprocess_lll(const matrix_of<B>& basis, LLLParams params = {})
{
    PreparedLLL<B> p;
    p.basis     = basis;
    p.transform = B::identity_int(B::cols(basis));
    p.gso       = gram_schmidt<B>(basis, params.gso);
    return p;
}

namespace detail {

// Lovász swap update (O(n) per swap, O(n²) total instead of O(n⁵) naïve).
//
// Swaps basis columns k-1 and k, then updates the GSO coefficient matrix and
// squared norms in-place without a full recompute.
//
// Derivation (Cohen, "A Course in Computational Algebraic Number Theory",
// Algorithm 2.6.3):
//
//   Let μ = mu[k][k-1], B = norms_sq[k-1], b = norms_sq[k], δ = b + μ²·B.
//
//   After the swap:
//     norms_sq[k-1]  ← δ
//     norms_sq[k]    ← B·b / δ
//     mu[k][k-1]     ← μ·B / δ
//
//   For rows i > k:
//     t1 = mu[i][k-1],  t2 = mu[i][k]
//     mu[i][k-1] ← (t2·b + t1·μ·B) / δ
//     mu[i][k]   ← t1 − μ·t2
//
//   For j < k-1: swap mu[k-1][j] and mu[k][j] (column re-index).
template <Backend B>
void lovasz_swap(matrix_of<B>& basis, int_matrix_of<B>& U,
                 GSOResult<B>& gso, std::size_t k)
{
    using S = scalar_of<B>;
    using I = integer_of<B>;

    const std::size_t n = B::cols(basis);
    const S mu_val = B::get(gso.mu, k, k - 1);
    const S B_val  = gso.norms_sq[k - 1];
    const S b_val  = gso.norms_sq[k];
    const S delta  = b_val + mu_val * mu_val * B_val;

    // Swap basis columns
    {
        vector_of<B> tmp = B::col(basis, k - 1);
        B::set_col(basis, k - 1, B::col(basis, k));
        B::set_col(basis, k, tmp);
    }

    // Swap unimodular transform columns (element-wise to avoid int_matrix column API)
    for (std::size_t i = 0; i < n; ++i) {
        const I tmp = B::int_get(U, i, k - 1);
        B::int_set(U, i, k - 1, B::int_get(U, i, k));
        B::int_set(U, i, k, tmp);
    }

    // Swap mu rows for j < k-1 (re-index the GSO columns)
    for (std::size_t j = 0; j < k - 1; ++j) {
        const S tmp = B::get(gso.mu, k - 1, j);
        B::set(gso.mu, k - 1, j, B::get(gso.mu, k, j));
        B::set(gso.mu, k, j, tmp);
    }

    // Lovász norm update
    gso.norms_sq[k - 1] = delta;
    gso.norms_sq[k]     = B_val * b_val / delta;
    B::set(gso.mu, k, k - 1, mu_val * B_val / delta);

    // Update mu rows above k
    for (std::size_t i = k + 1; i < n; ++i) {
        const S t1 = B::get(gso.mu, i, k - 1);
        const S t2 = B::get(gso.mu, i, k);
        B::set(gso.mu, i, k - 1, (t2 * b_val + t1 * mu_val * B_val) / delta);
        B::set(gso.mu, i, k,     t1 - mu_val * t2);
    }
}

} // namespace detail

// Stage 2: run LLL on a pre-computed PreparedLLL.
// Throws ConvergenceError if params.max_iter is exceeded.
template <Backend B = DefaultBackend>
[[nodiscard]] LLLResult<B>
lll_reduce(PreparedLLL<B> prepared, LLLParams params = {})
{
    using S = scalar_of<B>;
    using I = integer_of<B>;

    const std::size_t n = B::cols(prepared.basis);
    std::size_t k = 1;
    std::size_t swaps = 0;

    while (k < n) {
        // --- Size-reduction step ------------------------------------
        // For j from k-1 downto 0: if |mu(k,j)| > η, subtract a rounded
        // multiple of b_j from b_k and update mu and U accordingly.
        for (std::size_t ji = 0; ji < k; ++ji) {
            const std::size_t j = k - 1 - ji;
            const S mu_kj = B::get(prepared.gso.mu, k, j);

            if (std::abs(mu_kj) > static_cast<S>(params.eta)) {
                const I q  = static_cast<I>(std::round(mu_kj));
                const S qf = static_cast<S>(q);

                // b_k ← b_k − q·b_j
                B::set_col(prepared.basis, k,
                    B::subtract(B::col(prepared.basis, k),
                                B::scale(B::col(prepared.basis, j), qf)));

                // U_k ← U_k − q·U_j
                for (std::size_t i = 0; i < n; ++i)
                    B::int_set(prepared.transform, i, k,
                               B::int_get(prepared.transform, i, k)
                               - q * B::int_get(prepared.transform, i, j));

                // mu[k][l] ← mu[k][l] − q·mu[j][l]  for l < j
                for (std::size_t l = 0; l < j; ++l)
                    B::set(prepared.gso.mu, k, l,
                           B::get(prepared.gso.mu, k, l)
                           - qf * B::get(prepared.gso.mu, j, l));

                B::set(prepared.gso.mu, k, j, mu_kj - qf);
            }
        }

        // --- Lovász condition ----------------------------------------
        // ||b*_k||² ≥ (δ − mu(k,k-1)²)·||b*_{k-1}||²
        const S mu_k  = B::get(prepared.gso.mu, k, k - 1);
        const S delta = static_cast<S>(params.delta);
        const S lovasz_rhs = (delta - mu_k * mu_k) * prepared.gso.norms_sq[k - 1];

        if (prepared.gso.norms_sq[k] >= lovasz_rhs) {
            ++k;
        } else {
            if (params.max_iter > 0 && swaps >= params.max_iter)
                throw ConvergenceError("lll_reduce: swap limit exceeded");

            detail::lovasz_swap<B>(prepared.basis, prepared.transform, prepared.gso, k);
            ++swaps;

            if (k > 1) --k;
        }
    }

    LLLResult<B> result;
    result.reduced_basis = std::move(prepared.basis);
    result.transform     = std::move(prepared.transform);
    result.swap_count    = swaps;
    return result;
}

// Convenience overload: preprocess then reduce in a single call.
template <Backend B = DefaultBackend>
[[nodiscard]] LLLResult<B>
lll_reduce(const matrix_of<B>& basis, LLLParams params = {})
{
    return lll_reduce<B>(preprocess_lll<B>(basis, params), params);
}

} // namespace adrius
