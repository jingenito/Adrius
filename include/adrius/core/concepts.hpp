// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <concepts>
#include <cstddef>

namespace adrius {

// Minimum interface a linear algebra backend must satisfy.
//
// Design rules enforced here:
//   - Algorithms include only this header (and core/traits.hpp / core/params.hpp).
//   - No algorithm header may include any backend header directly.
//   - Basis vectors are COLUMNS of matrix_type (column-major convention).
//   - The unimodular transform acts on the right: M' = M · U.
//
// Adding a custom backend: implement a struct with all required type aliases
// and static member functions, then verify with static_assert(Backend<MyBackend>).
template <typename B>
concept Backend =
    requires {
        typename B::scalar_type;      // floating-point scalar (e.g. double)
        typename B::integer_type;     // signed integer        (e.g. int64_t)
        typename B::matrix_type;      // dense floating-point matrix
        typename B::vector_type;      // dense floating-point column vector
        typename B::int_matrix_type;  // dense integer matrix (unimodular transforms)
    }
    && std::floating_point<typename B::scalar_type>
    && std::signed_integral<typename B::integer_type>
    && requires(
        typename B::matrix_type&           M,
        const typename B::matrix_type&     cM,
        typename B::int_matrix_type&       iM,
        const typename B::int_matrix_type& ciM,
        const typename B::vector_type&     u,
        const typename B::vector_type&     v,
        typename B::scalar_type            s,
        typename B::integer_type           z,
        std::size_t r, std::size_t c, std::size_t k)
    {
        // Dimensions
        { B::rows(cM) } -> std::convertible_to<std::size_t>;
        { B::cols(cM) } -> std::convertible_to<std::size_t>;

        // Floating-point column access (basis vectors are columns)
        { B::col(cM, k) }       -> std::convertible_to<typename B::vector_type>;
        { B::set_col(M, k, u) } -> std::same_as<void>;

        // Floating-point element access (used for GSO coefficient matrix)
        { B::get(cM, r, c) }   -> std::convertible_to<typename B::scalar_type>;
        { B::set(M, r, c, s) } -> std::same_as<void>;

        // Integer element access (used for unimodular transform U)
        { B::int_get(ciM, r, c) }   -> std::convertible_to<typename B::integer_type>;
        { B::int_set(iM, r, c, z) } -> std::same_as<void>;

        // Vector arithmetic
        { B::inner_product(u, v) } -> std::convertible_to<typename B::scalar_type>;
        { B::norm_sq(v) }          -> std::convertible_to<typename B::scalar_type>;
        { B::scale(v, s) }         -> std::convertible_to<typename B::vector_type>;
        { B::add(u, v) }           -> std::convertible_to<typename B::vector_type>;
        { B::subtract(u, v) }      -> std::convertible_to<typename B::vector_type>;

        // Factory functions
        { B::make_zero_matrix(r, c) } -> std::convertible_to<typename B::matrix_type>;
        { B::make_zero_vector(r) }    -> std::convertible_to<typename B::vector_type>;
        { B::identity_int(r) }        -> std::convertible_to<typename B::int_matrix_type>;
    };

} // namespace adrius
