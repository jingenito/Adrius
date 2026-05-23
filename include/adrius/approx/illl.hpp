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
#include <span>
#include <vector>

namespace adrius {

// ── Result ────────────────────────────────────────────────────────────────

// Each element of `relations` is the integer vector (q, p₁, …, pₙ) found at
// one ILLL iteration.  `quality[k]` is max_i |q_k·αᵢ − p_{k,i}| for that
// iteration.  The Bosma-Smeets theorem guarantees the Dirichlet coefficient
// q^{1/n} · max|q·αᵢ − pᵢ| stays bounded by a constant depending only on n.
//
// Reference: W. Bosma & I. Smeets, "Finding simultaneous Diophantine
// approximations with prescribed quality", arXiv:1001.4455 (2010).
template <Backend B>
struct [[nodiscard]] ILLLResult {
    std::vector<std::vector<integer_of<B>>> relations;  // one (q,p) vector per iter
    std::vector<scalar_of<B>>              quality;     // max|q·αᵢ−pᵢ| per iter
    std::size_t                            iterations{0};

    ILLLResult() = default;
    ILLLResult(ILLLResult&&) noexcept = default;
    ILLLResult& operator=(ILLLResult&&) noexcept = default;
};

// ── Prepared state ────────────────────────────────────────────────────────

// Move-only lattice state for ILLL.  Holds the current (n+1)×(n+1) basis and
// the current scale c_k.  The basis is modified in-place between iterations.
template <Backend B>
struct PreparedILLL {
    matrix_of<B>              lattice;
    std::vector<scalar_of<B>> alpha;    // original targets (kept for extraction)
    scalar_of<B>              scale{};  // current c_k

    PreparedILLL() = default;
    PreparedILLL(PreparedILLL&&) noexcept = default;
    PreparedILLL& operator=(PreparedILLL&&) noexcept = default;
    PreparedILLL(const PreparedILLL&) = delete;
    PreparedILLL& operator=(const PreparedILLL&) = delete;
};

// ── Stage 1: preprocess ────────────────────────────────────────────────────
//
// Builds the Bosma-Smeets lattice for α = (α₁, …, αₙ):
//
//   col 0  = (c₀, α₁, α₂, …, αₙ)ᵀ   where  c₀ = ε^{n+1}
//   col j  = eⱼ  (identity block, j = 1…n)
//
// Key:  c₀ is chosen so that all columns have O(1) Euclidean norm, letting
// LLL balance the denominator direction against the error directions.
// Contrast with the naïve construction (N, N·α₁, …, Nαₙ) which has
// c₀-column norm ≈ N ≫ 1, causing LLL to ignore it entirely (q = 0).
template <Backend B = DefaultBackend>
[[nodiscard]] PreparedILLL<B>
preprocess_illl(std::span<const scalar_of<B>> alpha, ILLLParams params = {})
{
    if (alpha.empty())
        throw DomainError("preprocess_illl: alpha must be non-empty");
    if (params.epsilon <= 0.0 || params.epsilon >= 1.0)
        throw DomainError("preprocess_illl: epsilon must be in (0, 1)");

    const std::size_t n   = alpha.size();
    const std::size_t dim = n + 1;

    // Initial scale: c₀ = ε^{n+1}  (Bosma-Smeets Algorithm 1, k = 0)
    const scalar_of<B> c0 = static_cast<scalar_of<B>>(
        std::pow(params.epsilon, static_cast<double>(n + 1)));

    PreparedILLL<B> p;
    p.alpha.assign(alpha.begin(), alpha.end());
    p.scale   = c0;
    p.lattice = B::make_zero_matrix(dim, dim);

    // Column 0: the denominator direction scaled by c₀, targets as-is
    B::set(p.lattice, 0, 0, c0);
    for (std::size_t i = 0; i < n; ++i)
        B::set(p.lattice, i + 1, 0, alpha[i]);

    // Columns 1..n: identity block
    for (std::size_t j = 1; j < dim; ++j)
        B::set(p.lattice, j, j, scalar_of<B>{1});

    return p;
}

// ── Stage 2: iterate ───────────────────────────────────────────────────────
//
// Outer loop (Bosma & Smeets 2010, Algorithm 1):
//
//   1. LLL-reduce the current lattice.
//   2. The first column b of the reduced basis encodes the approximation:
//        b[0]   = q · c_k            (denominator, scaled)
//        b[i>0] = q · αᵢ − pᵢ       (approximation errors, directly)
//   3. Extract q and pᵢ; record quality = max_i |b[i>0]|.
//   4. Scale update (Lemma 3.4 / Theorem 3.5, eq. 24-28):
//        c_{k+1} = c_k · 2^{−n(n+1)/4}
//   5. Rescale column 0 of the current reduced basis by c_{k+1}/c_k.
//      The resulting basis is nearly reduced for the new lattice, so
//      the next LLL call is efficient.
//
// Terminates when q ≥ max_denominator, quality < quality_tol, or
// max_iterations is reached.
template <Backend B = DefaultBackend>
[[nodiscard]] ILLLResult<B>
illl(PreparedILLL<B> prepared, ILLLParams params = {})
{
    // ADL: for built-in types finds std::abs/round; for Boost MP types finds
    // boost::multiprecision::abs/round via argument-dependent lookup.
    using std::abs;
    using std::round;

    const std::size_t n   = prepared.alpha.size();
    const std::size_t dim = n + 1;

    // Per-step scale factor: c(k+1) = c(k) · 2^{-n(n+1)/4}
    // Derivation: from c(k) = (2^{-kn/4}·ε)^{n+1}, the ratio of consecutive
    // terms is 2^{-n/4·(n+1)} = 2^{-n(n+1)/4}.
    const scalar_of<B> scale_ratio = static_cast<scalar_of<B>>(
        std::pow(2.0, -static_cast<double>(n * (n + 1)) / 4.0));

    ILLLResult<B> result;
    result.relations.reserve(params.max_iterations);
    result.quality.reserve(params.max_iterations);

    for (std::size_t iter = 0; iter < params.max_iterations; ++iter) {
        // Bail if scale has dropped so far that gram_schmidt will fail.
        // The basis becomes ill-conditioned and column 0 is essentially zero.
        if (iter > 0 && prepared.scale < 1e-14)
            break;

        auto lll = lll_reduce<B>(prepared.lattice, params.lll);

        // b[0] = q · c_k  →  q = round(|b[0]| / c_k)
        const scalar_of<B> b0 = B::get(lll.reduced_basis, 0, 0);
        const integer_of<B> q =
            static_cast<integer_of<B>>(round(abs(b0) / prepared.scale));

        if (q > integer_of<B>{0}) {
            // b[i] = q·αᵢ − pᵢ  →  pᵢ = round(q·αᵢ − b[i])
            std::vector<integer_of<B>> relation(dim);
            relation[0] = q;
            scalar_of<B> max_err{0};
            for (std::size_t i = 0; i < n; ++i) {
                const scalar_of<B> bi =
                    B::get(lll.reduced_basis, i + 1, 0);
                relation[i + 1] = static_cast<integer_of<B>>(
                    round(static_cast<scalar_of<B>>(q) * prepared.alpha[i] - bi));
                scalar_of<B> abs_bi = abs(bi);
                if (abs_bi > max_err) max_err = abs_bi;
            }
            result.relations.push_back(std::move(relation));
            result.quality.push_back(max_err);

            if (params.max_denominator > 0 &&
                q >= static_cast<integer_of<B>>(params.max_denominator))
                break;
            if (max_err < static_cast<scalar_of<B>>(params.quality_tol))
                break;
        }

        // Rescale column 0 of the LLL-reduced basis by c_{k+1}/c_k.
        // The other columns are unchanged.  This leaves the basis nearly
        // reduced for the new lattice, amortising the cost of the next LLL.
        prepared.lattice = std::move(lll.reduced_basis);
        for (std::size_t i = 0; i < dim; ++i)
            B::set(prepared.lattice, i, 0,
                   B::get(prepared.lattice, i, 0) * scale_ratio);
        prepared.scale *= scale_ratio;
    }

    result.iterations = result.relations.size();
    return result;
}

// ── Convenience overload ─────────────────────────────────────────────────
template <Backend B = DefaultBackend>
[[nodiscard]] ILLLResult<B>
illl(std::span<const scalar_of<B>> alpha, ILLLParams params = {})
{
    return illl<B>(preprocess_illl<B>(alpha, params), params);
}

} // namespace adrius
