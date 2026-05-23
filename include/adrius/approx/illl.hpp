// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/core/concepts.hpp>
#include <adrius/core/error.hpp>
#include <adrius/core/params.hpp>
#include <adrius/core/traits.hpp>
#include <adrius/backend/default_backend.hpp>
#include <adrius/linalg/lll.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace adrius {

// Result of ILLL (Iterated LLL) from Bosma & Smeets (2010).
//
//   relations:  integer matrix whose columns are the simultaneous approximation
//               vectors found; column j gives (q_j, p_{j,1}, …, p_{j,n})
//   quality:    ||q·α − p||_∞ · q for each relation (smaller is better)
//   iterations: number of outer ILLL iterations performed
//
// Reference: W. Bosma & I. Smeets, "Finding simultaneous Diophantine
// approximations with prescribed quality", 2010.
template <Backend B>
struct [[nodiscard]] ILLLResult {
    int_matrix_of<B>           relations;
    std::vector<scalar_of<B>>  quality;
    std::size_t                iterations{0};

    ILLLResult() = default;
    ILLLResult(ILLLResult&&) noexcept = default;
    ILLLResult& operator=(ILLLResult&&) noexcept = default;
};

// Pre-computed state for ILLL. Move-only — lattice matrices are large.
//
// ILLL lattice construction for α = (α_1, …, α_n):
// Build the (n+1)×(n+1) matrix at an initial scale N_0:
//
//   [  N_0    0      0    … 0  ]
//   [ N_0·α_1 1      0    … 0  ]
//   [ N_0·α_2 0      1    … 0  ]
//   [   ⋮                  ⋮   ]
//   [ N_0·α_n 0      0    … 1  ]
//
// After each LLL reduction the scale N is updated and the first column is
// rescaled, implementing the "iterated" part of ILLL.
template <Backend B>
struct PreparedILLL {
    matrix_of<B>              lattice;    // current working lattice
    int_matrix_of<B>          transform;  // cumulative unimodular transform
    std::vector<scalar_of<B>> alpha;      // approximation target
    scalar_of<B>              scale{};    // current N

    PreparedILLL() = default;
    PreparedILLL(PreparedILLL&&) noexcept = default;
    PreparedILLL& operator=(PreparedILLL&&) noexcept = default;
    PreparedILLL(const PreparedILLL&) = delete;
    PreparedILLL& operator=(const PreparedILLL&) = delete;
};

// Stage 1: construct the initial ILLL lattice for α at scale N_0.
// alpha must be non-empty and N_0 must be positive.
template <Backend B = DefaultBackend>
[[nodiscard]] PreparedILLL<B>
preprocess_illl(std::span<const scalar_of<B>> alpha,
                scalar_of<B> initial_scale = scalar_of<B>{1e6})
{
    if (alpha.empty())
        throw DomainError("preprocess_illl: alpha must be non-empty");
    if (initial_scale <= scalar_of<B>{0})
        throw DomainError("preprocess_illl: initial_scale must be positive");

    const std::size_t n   = alpha.size();
    const std::size_t dim = n + 1;

    PreparedILLL<B> p;
    p.alpha.assign(alpha.begin(), alpha.end());
    p.scale   = initial_scale;
    p.lattice = B::make_zero_matrix(dim, dim);

    // First column: (N_0, ⌊N_0·α_1⌋, …, ⌊N_0·α_n⌋)ᵀ
    B::set(p.lattice, 0, 0, initial_scale);
    for (std::size_t i = 0; i < n; ++i)
        B::set(p.lattice, i + 1, 0,
               static_cast<scalar_of<B>>(std::round(initial_scale * alpha[i])));

    // Identity block for the remaining columns
    for (std::size_t j = 1; j < dim; ++j)
        B::set(p.lattice, j, j, scalar_of<B>{1});

    p.transform = B::identity_int(dim);
    return p;
}

// Stage 2: run ILLL on the pre-built PreparedILLL.
//
// The outer loop (Bosma & Smeets Algorithm 1):
//   1. LLL-reduce the current lattice.
//   2. Extract the shortest vector v_0; compute quality Q_0 = ||q·α − p||_∞·q.
//   3. Compute the new scale N_{k+1} from the current reduction quality
//      (implements eqs. 24 and 28 from Bosma & Smeets 2010; see decision record
//      §Key Algorithms for the bounds that the previous MUSE implementation
//      got wrong).
//   4. Rescale the first lattice column by N_{k+1}/N_k and iterate.
//
// TODO: Implement the scale update rule from Theorem 23 (eq. 24) and
//       Lemma 25 (eq. 28) of Bosma & Smeets (2010). These bounds must be
//       derived carefully — do not estimate from the previous implementation.
template <Backend B = DefaultBackend>
[[nodiscard]] ILLLResult<B>
illl(PreparedILLL<B> prepared, ILLLParams params = {})
{
    const std::size_t n   = prepared.alpha.size();
    const std::size_t dim = n + 1;

    ILLLResult<B> result;
    result.relations = B::identity_int(dim);
    result.quality.reserve(params.max_iterations > 0 ? params.max_iterations : 16);

    scalar_of<B> prev_quality = std::numeric_limits<scalar_of<B>>::max();
    std::size_t  iter         = 0;
    const std::size_t limit   = params.max_iterations > 0
                                    ? params.max_iterations
                                    : std::size_t{1000};

    while (iter < limit) {
        // LLL-reduce the current lattice
        auto lll = lll_reduce<B>(prepared.lattice, params.lll);

        // Extract approximation from the first column of the reduced basis
        const scalar_of<B> q_scaled = B::get(lll.reduced_basis, 0, 0);
        const integer_of<B> q =
            static_cast<integer_of<B>>(std::round(q_scaled / prepared.scale));

        scalar_of<B> max_err{0};
        for (std::size_t i = 0; i < n; ++i) {
            const integer_of<B> p =
                static_cast<integer_of<B>>(std::round(B::get(lll.reduced_basis, i + 1, 0)));
            const scalar_of<B> err =
                std::abs(static_cast<scalar_of<B>>(q) * prepared.alpha[i]
                         - static_cast<scalar_of<B>>(p));
            max_err = std::max(max_err, err);
        }
        const scalar_of<B> quality = max_err * static_cast<scalar_of<B>>(std::abs(q));
        result.quality.push_back(quality);

        // Convergence check
        const scalar_of<B> improvement =
            (prev_quality - quality) / (prev_quality + scalar_of<B>{1e-300});
        if (improvement < static_cast<scalar_of<B>>(params.convergence_tol) && iter > 0)
            break;

        prev_quality = quality;

        // TODO: implement the Bosma-Smeets scale update (Theorem 23, eq. 24;
        //       Lemma 25, eq. 28). The new scale N_{k+1} is derived from the
        //       norms of the LLL-reduced basis. Placeholder: geometric growth.
        const scalar_of<B> new_scale = prepared.scale * static_cast<scalar_of<B>>(10);

        // Rescale the first lattice column by N_{k+1} / N_k
        const scalar_of<B> ratio = new_scale / prepared.scale;
        prepared.lattice = std::move(lll.reduced_basis);
        for (std::size_t i = 0; i < dim; ++i)
            B::set(prepared.lattice, i, 0, B::get(prepared.lattice, i, 0) * ratio);
        prepared.scale = new_scale;

        ++iter;
    }

    result.iterations = iter;
    return result;
}

// Convenience overload.
template <Backend B = DefaultBackend>
[[nodiscard]] ILLLResult<B>
illl(std::span<const scalar_of<B>> alpha,
     ILLLParams params = {},
     scalar_of<B> initial_scale = scalar_of<B>{1e6})
{
    return illl<B>(preprocess_illl<B>(alpha, initial_scale), params);
}

} // namespace adrius
