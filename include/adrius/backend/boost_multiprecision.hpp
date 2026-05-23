// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/core/concepts.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace adrius {

// ── Boost.Multiprecision Backend ──────────────────────────────────────────
//
// Provides arbitrary-precision arithmetic via Boost.Multiprecision.
// Useful when Eigen's double precision (15-17 significant digits) is insufficient
// for ill-conditioned problems or high-denominator approximations.
//
// Usage:
//   // 100 decimal digits of precision (default: 50)
//   using HighPrecision = adrius::BoostBackend<100>;
//   auto result = adrius::illl<HighPrecision>(alpha, params);
//
// Type Parameters:
//   Digits: Decimal digits of precision for floating-point (default 50)
//   - 50 ≈ 166 bits (enough for most uses)
//   - 100 ≈ 332 bits (high-precision computations)
//   - 1000+ for extreme cases (slower)

template <unsigned Digits = 50>
struct BoostBackend {
    // Floating-point type: decimal float with configurable precision
    using scalar_type = boost::multiprecision::number<
        boost::multiprecision::cpp_dec_float<Digits>>;

    // Integer type: arbitrary-precision signed integer
    using integer_type = boost::multiprecision::cpp_int;

    // Dense matrix over scalar_type
    using matrix_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;

    // Dense vector
    using vector_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

    // Integer matrix for unimodular transforms
    using int_matrix_type = Eigen::Matrix<integer_type, Eigen::Dynamic, Eigen::Dynamic>;

    // ── Basic matrix operations ───────────────────────────────────────────

    static std::size_t rows(const matrix_type& m) { return m.rows(); }
    static std::size_t cols(const matrix_type& m) { return m.cols(); }

    static scalar_type get(const matrix_type& m, std::size_t i, std::size_t j) {
        return m(i, j);
    }

    static void set(matrix_type& m, std::size_t i, std::size_t j, const scalar_type& val) {
        m(i, j) = val;
    }

    // ── Integer matrix operations ─────────────────────────────────────────

    static integer_type get_int(const int_matrix_type& m, std::size_t i, std::size_t j) {
        return m(i, j);
    }

    static void set_int(int_matrix_type& m, std::size_t i, std::size_t j,
                        const integer_type& val) {
        m(i, j) = val;
    }

    // ── Factory functions ─────────────────────────────────────────────────

    static matrix_type make_zero_matrix(std::size_t rows, std::size_t cols) {
        return matrix_type::Zero(rows, cols);
    }

    static vector_type make_zero_vector(std::size_t size) {
        return vector_type::Zero(size);
    }

    static int_matrix_type make_zero_int_matrix(std::size_t rows, std::size_t cols) {
        return int_matrix_type::Zero(rows, cols);
    }

    static int_matrix_type make_identity_int_matrix(std::size_t size) {
        return int_matrix_type::Identity(size, size);
    }

    // Map a std::span into an Eigen vector (zero-copy)
    static vector_type map_vector(std::span<const scalar_type> data) {
        if (data.empty()) throw std::runtime_error("Cannot map empty span");
        return Eigen::Map<const vector_type>(data.data(), data.size());
    }

    // ── Vector operations ─────────────────────────────────────────────────

    static scalar_type inner_product(const vector_type& a, const vector_type& b) {
        if (a.size() != b.size())
            throw std::runtime_error(
                "inner_product: size mismatch " + std::to_string(a.size()) + " vs " +
                std::to_string(b.size()));
        return a.dot(b);
    }

    static scalar_type squared_norm(const vector_type& v) { return v.squaredNorm(); }

    static scalar_type norm(const vector_type& v) { return v.norm(); }

    static vector_type scale(const vector_type& v, const scalar_type& s) { return v * s; }

    static vector_type add(const vector_type& a, const vector_type& b) { return a + b; }

    static vector_type subtract(const vector_type& a, const vector_type& b) {
        return a - b;
    }

    // ── Conversions ───────────────────────────────────────────────────────

    // Convert scalar to double (may lose precision if value uses full width)
    static double to_double(const scalar_type& val) {
        return static_cast<double>(val);
    }

    // Convert scalar to integer (rounds toward zero)
    static integer_type to_integer(const scalar_type& val) {
        return integer_type(val);  // Truncates toward zero
    }

    // Convert integer to scalar
    static scalar_type from_integer(const integer_type& val) { return scalar_type(val); }

    // ── Utility ───────────────────────────────────────────────────────────

    // Absolute value
    static scalar_type abs(const scalar_type& val) { return boost::multiprecision::abs(val); }

    // Maximum of two values
    static scalar_type max(const scalar_type& a, const scalar_type& b) {
        return a > b ? a : b;
    }

    // Minimum of two values
    static scalar_type min(const scalar_type& a, const scalar_type& b) {
        return a < b ? a : b;
    }

    // Square root
    static scalar_type sqrt(const scalar_type& val) {
        if (val < 0)
            throw std::domain_error("sqrt of negative number");
        return boost::multiprecision::sqrt(val);
    }

    // Floor function
    static scalar_type floor(const scalar_type& val) {
        return boost::multiprecision::floor(val);
    }

    // Ceiling function
    static scalar_type ceil(const scalar_type& val) {
        return boost::multiprecision::ceil(val);
    }

    // Round to nearest integer
    static scalar_type round(const scalar_type& val) {
        return boost::multiprecision::round(val);
    }

    // Power function
    static scalar_type pow(const scalar_type& base, double exponent) {
        return boost::multiprecision::pow(base, exponent);
    }
};

// ── Type aliases for common precisions ────────────────────────────────────

// Standard precision: 50 decimal digits (166 bits, ~1e-50 error)
using BoostBackendDefault = BoostBackend<50>;

// High precision: 100 decimal digits (332 bits, ~1e-100 error)
using BoostBackendHighPrecision = BoostBackend<100>;

// Very high precision: 200 decimal digits (664 bits, ~1e-200 error)
using BoostBackendVeryHighPrecision = BoostBackend<200>;

// ── Static assertion for concept compliance ───────────────────────────────

static_assert(Backend<BoostBackendDefault>,
              "BoostBackendDefault must satisfy Backend concept");

}  // namespace adrius
