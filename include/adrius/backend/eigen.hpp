// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

// ============================================================
// This is the ONLY header in the library permitted to include
// Eigen headers. All algorithms are written against the Backend
// concept in core/concepts.hpp; they must never include this
// file directly.
// ============================================================

#include <adrius/core/concepts.hpp>

// Suppress warnings from Eigen's own headers. This is the one place in the
// library where platform-specific pragmas are unavoidable — they guard a
// third-party include, not our own code.
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
#include <cstdint>
#include <span>

namespace adrius {

struct EigenBackend {
    using scalar_type     = double;
    using integer_type    = std::int64_t;
    using matrix_type     = Eigen::MatrixXd;
    using vector_type     = Eigen::VectorXd;
    using int_matrix_type = Eigen::Matrix<std::int64_t, Eigen::Dynamic, Eigen::Dynamic>;

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
                                          std::size_t r, std::size_t c) noexcept {
        return M(idx(r), idx(c));
    }
    static void set(matrix_type& M, std::size_t r, std::size_t c, scalar_type s) noexcept {
        M(idx(r), idx(c)) = s;
    }

    // --- Integer element access (unimodular transform U) ------------
    [[nodiscard]] static integer_type int_get(const int_matrix_type& M,
                                               std::size_t r, std::size_t c) noexcept {
        return M(idx(r), idx(c));
    }
    static void int_set(int_matrix_type& M,
                        std::size_t r, std::size_t c, integer_type z) noexcept {
        M(idx(r), idx(c)) = z;
    }

    // --- Vector arithmetic ------------------------------------------
    [[nodiscard]] static scalar_type inner_product(const vector_type& u,
                                                    const vector_type& v) noexcept {
        return u.dot(v);
    }
    [[nodiscard]] static scalar_type norm_sq(const vector_type& v) noexcept {
        return v.squaredNorm();
    }
    [[nodiscard]] static vector_type scale(const vector_type& v, scalar_type s) {
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

    // --- Zero-copy views of user-provided memory (column-major) -----
    // The returned Map is valid only while the source buffer is alive.
    [[nodiscard]] static Eigen::Map<const vector_type>
    map_vector(const scalar_type* data, std::size_t n) noexcept {
        return {data, idx(n)};
    }
    [[nodiscard]] static Eigen::Map<const vector_type>
    map_vector(std::span<const scalar_type> s) noexcept {
        return {s.data(), idx(s.size())};
    }
    [[nodiscard]] static Eigen::Map<const matrix_type>
    map_matrix(const scalar_type* data, std::size_t r, std::size_t c) noexcept {
        return {data, idx(r), idx(c)};
    }

private:
    static constexpr Eigen::Index idx(std::size_t n) noexcept {
        return static_cast<Eigen::Index>(n);
    }
};

static_assert(Backend<EigenBackend>,
    "EigenBackend must satisfy the Backend concept — a required member is missing");

} // namespace adrius
