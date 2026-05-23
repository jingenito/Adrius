// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

// ============================================================
// This is the ONLY header permitted to include Boost headers.
// Algorithms must never include this file directly — they use
// only core/concepts.hpp and the Backend abstraction.
// ============================================================

#include <adrius/core/concepts.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#if defined(_MSC_VER)
#  pragma warning(push, 0)
#elif defined(__GNUC__) || defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wall"
#  pragma GCC diagnostic ignored "-Wextra"
#endif

#include <Eigen/Dense>

#if defined(_MSC_VER)
#  pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#  pragma GCC diagnostic pop
#endif
#include <cstddef>
#include <span>

namespace adrius {

// ── Boost.Multiprecision Backend ──────────────────────────────────────────
//
// Arbitrary-precision arithmetic via Boost.Multiprecision. Use when Eigen's
// double precision (15-17 digits) is insufficient for ill-conditioned problems
// or high-denominator approximations.
//
// Usage:
//   using HP = adrius::BoostBackend<100>;   // 100 decimal digits
//   auto result = adrius::illl<HP>(alpha, params);
//
// Template parameter:
//   Digits — decimal digits of precision for the floating-point scalar type.
//             50 ≈ 166 bits, 100 ≈ 332 bits.

template <unsigned Digits = 50>
struct BoostBackend {
    using scalar_type     = boost::multiprecision::number<
                                boost::multiprecision::cpp_dec_float<Digits>>;
    using integer_type    = boost::multiprecision::cpp_int;
    using matrix_type     = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
    using vector_type     = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using int_matrix_type = Eigen::Matrix<integer_type, Eigen::Dynamic, Eigen::Dynamic>;

    // --- Dimensions -------------------------------------------------
    [[nodiscard]] static std::size_t rows(const matrix_type& M) noexcept {
        return static_cast<std::size_t>(M.rows());
    }
    [[nodiscard]] static std::size_t cols(const matrix_type& M) noexcept {
        return static_cast<std::size_t>(M.cols());
    }

    // --- Column access ----------------------------------------------
    [[nodiscard]] static vector_type col(const matrix_type& M, std::size_t k) {
        return M.col(idx(k));
    }
    static void set_col(matrix_type& M, std::size_t k, const vector_type& v) {
        M.col(idx(k)) = v;
    }

    // --- Floating-point element access (GSO coefficient matrix) -----
    [[nodiscard]] static scalar_type get(const matrix_type& M,
                                          std::size_t r, std::size_t c) {
        return M(idx(r), idx(c));
    }
    static void set(matrix_type& M, std::size_t r, std::size_t c, const scalar_type& s) {
        M(idx(r), idx(c)) = s;
    }

    // --- Integer element access (unimodular transform U) ------------
    [[nodiscard]] static integer_type int_get(const int_matrix_type& M,
                                               std::size_t r, std::size_t c) {
        return M(idx(r), idx(c));
    }
    static void int_set(int_matrix_type& M,
                        std::size_t r, std::size_t c, const integer_type& z) {
        M(idx(r), idx(c)) = z;
    }

    // --- Vector arithmetic ------------------------------------------
    [[nodiscard]] static scalar_type inner_product(const vector_type& u,
                                                    const vector_type& v) {
        return u.dot(v);
    }
    [[nodiscard]] static scalar_type norm_sq(const vector_type& v) {
        return v.squaredNorm();
    }
    [[nodiscard]] static vector_type scale(const vector_type& v, const scalar_type& s) {
        return s * v;
    }
    [[nodiscard]] static vector_type add(const vector_type& u, const vector_type& v) {
        return u + v;
    }
    [[nodiscard]] static vector_type subtract(const vector_type& u, const vector_type& v) {
        return u - v;
    }

    // --- Factory functions ------------------------------------------
    [[nodiscard]] static matrix_type make_zero_matrix(std::size_t r, std::size_t c) {
        return matrix_type::Zero(idx(r), idx(c));
    }
    [[nodiscard]] static vector_type make_zero_vector(std::size_t r) {
        return vector_type::Zero(idx(r));
    }
    [[nodiscard]] static int_matrix_type identity_int(std::size_t r) {
        return int_matrix_type::Identity(idx(r), idx(r));
    }

    // --- Zero-copy view of user-provided data -----------------------
    [[nodiscard]] static Eigen::Map<const vector_type>
    map_vector(std::span<const scalar_type> s) noexcept {
        return {s.data(), idx(s.size())};
    }

private:
    static constexpr Eigen::Index idx(std::size_t n) noexcept {
        return static_cast<Eigen::Index>(n);
    }
};

// ── Type aliases for common precisions ────────────────────────────────────

using BoostBackendDefault         = BoostBackend<50>;   // ~1e-50 error
using BoostBackendHighPrecision   = BoostBackend<100>;  // ~1e-100 error
using BoostBackendVeryHighPrecision = BoostBackend<200>; // ~1e-200 error

// ── Static assertion for concept compliance ───────────────────────────────

static_assert(Backend<BoostBackendDefault>,
              "BoostBackendDefault must satisfy the Backend concept");

} // namespace adrius
