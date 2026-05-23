// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/core/error.hpp>

#include <compare>
#include <concepts>
#include <cstdint>
#include <numeric>
#include <ostream>

namespace adrius {

// Exact rational number stored in lowest terms with denominator > 0.
// Int must be a signed integer type; defaults to int64_t.
//
// All arithmetic operators return a new reduced Rational — no mutation.
// constexpr throughout so compile-time rational arithmetic is possible.
template <std::signed_integral Int = std::int64_t>
class Rational {
public:
    constexpr Rational() noexcept : num_(0), den_(1) {}
    constexpr explicit Rational(Int n) noexcept : num_(n), den_(1) {}
    constexpr Rational(Int numerator, Int denominator) : num_(numerator), den_(denominator) {
        normalize();
    }

    [[nodiscard]] constexpr Int numerator()   const noexcept { return num_; }
    [[nodiscard]] constexpr Int denominator() const noexcept { return den_; }

    [[nodiscard]] constexpr double to_double() const noexcept {
        return static_cast<double>(num_) / static_cast<double>(den_);
    }

    [[nodiscard]] constexpr bool is_integer() const noexcept { return den_ == Int{1}; }
    [[nodiscard]] constexpr bool is_zero()    const noexcept { return num_ == Int{0}; }

    // Arithmetic — results are always reduced.
    [[nodiscard]] constexpr Rational operator+(const Rational& rhs) const {
        return {num_ * rhs.den_ + rhs.num_ * den_, den_ * rhs.den_};
    }
    [[nodiscard]] constexpr Rational operator-(const Rational& rhs) const {
        return {num_ * rhs.den_ - rhs.num_ * den_, den_ * rhs.den_};
    }
    [[nodiscard]] constexpr Rational operator*(const Rational& rhs) const {
        return {num_ * rhs.num_, den_ * rhs.den_};
    }
    [[nodiscard]] constexpr Rational operator/(const Rational& rhs) const {
        return {num_ * rhs.den_, den_ * rhs.num_};
    }
    [[nodiscard]] constexpr Rational operator-() const noexcept {
        Rational r; r.num_ = -num_; r.den_ = den_; return r;
    }

    constexpr Rational& operator+=(const Rational& rhs) { return *this = *this + rhs; }
    constexpr Rational& operator-=(const Rational& rhs) { return *this = *this - rhs; }
    constexpr Rational& operator*=(const Rational& rhs) { return *this = *this * rhs; }
    constexpr Rational& operator/=(const Rational& rhs) { return *this = *this / rhs; }

    // Comparison. Cross-multiplication avoids floating-point; both denominators
    // are positive by invariant so the sign of the product is correct.
    [[nodiscard]] constexpr auto operator<=>(const Rational& rhs) const noexcept {
        return (num_ * rhs.den_) <=> (rhs.num_ * den_);
    }
    [[nodiscard]] constexpr bool operator==(const Rational&) const noexcept = default;

    friend std::ostream& operator<<(std::ostream& os, const Rational& r) {
        if (r.den_ == Int{1}) return os << r.num_;
        return os << r.num_ << '/' << r.den_;
    }

private:
    Int num_;
    Int den_;

    constexpr void normalize() {
        if (den_ == Int{0})
            throw DomainError("Rational: zero denominator");
        if (den_ < Int{0}) { num_ = -num_; den_ = -den_; }
        const Int g = std::gcd(num_ < Int{0} ? -num_ : num_, den_);
        if (g > Int{1}) { num_ /= g; den_ /= g; }
    }
};

} // namespace adrius
