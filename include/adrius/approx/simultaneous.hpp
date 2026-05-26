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

// Result of simultaneous Diophantine approximation:
// finds q ∈ ℤ and p ∈ ℤⁿ such that |q·αᵢ − pᵢ| is small for all i.
//
//   denominator: the common denominator q found
//   numerators:  the integer vector p = (p_1, …, p_n)
//   quality:     max_i |q·αᵢ − pᵢ|  (smaller is better)
template <Backend B>
struct [[nodiscard]] SimultApproxResult {
    integer_of<B>              denominator{};
    std::vector<integer_of<B>> numerators;
    scalar_of<B>               quality{};

    SimultApproxResult() = default;
    SimultApproxResult(SimultApproxResult&&) noexcept = default;
    SimultApproxResult& operator=(SimultApproxResult&&) noexcept = default;
};

// Pre-computed state for simultaneous approximation.
// Holds the constructed lattice for α at scale N. Move-only.
//
// Lattice construction (standard same-denominator approach):
// Given α = (α_1, …, α_n) and scale N, build the (n+1)×(n+1) matrix:
//
//   [ N    0   0  … 0 ]
//   [ Nα_1 1   0  … 0 ]
//   [ Nα_2 0   1  … 0 ]
//   [  ⋮        ⋱    ⋮ ]
//   [ Nα_n 0   0  … 1 ]
//
// The first column of the LLL-reduced basis yields (q, p_1, …, p_n) with
// |q·αᵢ − pᵢ| ≈ 1/N^{1/n} by Dirichlet's theorem.
template <Backend B>
struct PreparedSimultaneous {
    matrix_of<B>              lattice;
    std::vector<scalar_of<B>> alpha;   // original target (for quality computation)
    scalar_of<B>              scale{}; // N

    PreparedSimultaneous() = default;
    PreparedSimultaneous(PreparedSimultaneous&&) noexcept = default;
    PreparedSimultaneous& operator=(PreparedSimultaneous&&) noexcept = default;
    PreparedSimultaneous(const PreparedSimultaneous&) = delete;
    PreparedSimultaneous& operator=(const PreparedSimultaneous&) = delete;
};

// Stage 1: build the approximation lattice for α at scale N.
// alpha must be non-empty; N must be positive.
template <Backend B = DefaultBackend>
[[nodiscard]] PreparedSimultaneous<B>
preprocess_simultaneous(std::span<const scalar_of<B>> alpha,
                        scalar_of<B> N)
{
    // ADL: finds boost::multiprecision::round for MP types, std:: for built-ins.
    using std::round;

    if (alpha.empty())
        throw DomainError("preprocess_simultaneous: alpha must be non-empty");
    if (N <= scalar_of<B>{0})
        throw DomainError("preprocess_simultaneous: scale N must be positive");

    const std::size_t n = alpha.size();
    const std::size_t dim = n + 1;

    PreparedSimultaneous<B> p;
    p.alpha.assign(alpha.begin(), alpha.end());
    p.scale  = N;
    p.lattice = B::make_zero_matrix(dim, dim);

    // First column: (N, N·α_1, …, N·α_n)ᵀ
    B::set(p.lattice, 0, 0, N);
    for (std::size_t i = 0; i < n; ++i)
        B::set(p.lattice, i + 1, 0, round(N * alpha[i]));

    // Remaining columns: identity block
    for (std::size_t j = 1; j < dim; ++j)
        B::set(p.lattice, j, j, scalar_of<B>{1});

    return p;
}

// Stage 2: LLL-reduce the pre-built lattice and extract the approximation.
template <Backend B = DefaultBackend>
[[nodiscard]] SimultApproxResult<B>
simultaneous_approx(PreparedSimultaneous<B> prepared, LLLParams params = {})
{
    // ADL: finds boost::multiprecision::abs/round for MP types, std:: for built-ins.
    using std::abs;
    using std::round;

    auto lll_result = lll_reduce<B>(prepared.lattice, params);
    const std::size_t n = prepared.alpha.size();

    // The first column of the reduced basis: (q_scaled, p_1, …, p_n)
    const scalar_of<B> q_scaled = B::get(lll_result.reduced_basis, 0, 0);
    const integer_of<B> q = static_cast<integer_of<B>>(round(q_scaled / prepared.scale));

    SimultApproxResult<B> result;
    result.denominator = q;
    result.numerators.resize(n);

    scalar_of<B> max_err{0};
    for (std::size_t i = 0; i < n; ++i) {
        result.numerators[i] =
            static_cast<integer_of<B>>(round(B::get(lll_result.reduced_basis, i + 1, 0)));
        scalar_of<B> err = abs(
            static_cast<scalar_of<B>>(q) * prepared.alpha[i]
            - static_cast<scalar_of<B>>(result.numerators[i]));
        if (err > max_err) max_err = err;
    }
    result.quality = max_err;

    return result;
}

// Convenience overload.
template <Backend B = DefaultBackend>
[[nodiscard]] SimultApproxResult<B>
simultaneous_approx(std::span<const scalar_of<B>> alpha,
                    scalar_of<B> N,
                    LLLParams params = {})
{
    return simultaneous_approx<B>(preprocess_simultaneous<B>(alpha, N), params);
}

} // namespace adrius
