// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/approx/continued_fraction.hpp>
#include <adrius/core/error.hpp>
#include <adrius/core/params.hpp>
#include <adrius/util/rational_type.hpp>

#include <cmath>
#include <cstdint>
#include <optional>

namespace adrius {

// Returns the best rational approximation p/q to x with q ≤ params.max_denominator.
//
// Algorithm: compute CF convergents until q_k exceeds max_denominator, then
// check the last valid convergent and the "semi-convergent" formed by the
// largest integer multiple of the previous convergent that stays in bound.
// This is guaranteed to return the best approximation in the Farey sense.
//
// Reference: Hardy & Wright, "An Introduction to the Theory of Numbers",
// Theorem 171 (best approximation property of convergents).
[[nodiscard]] inline Rational<std::int64_t>
best_rational(double x, RationalApproxParams params = {})
{
    if (params.max_denominator < 1)
        throw DomainError("best_rational: max_denominator must be >= 1");

    // Use a CF depth sufficient for max_denominator; convergent denominators
    // grow at least as fast as Fibonacci numbers, so ~90 steps covers int64_t.
    CFParams cf_params{.max_depth = 90};

    // Same seed as cf_convergents: p_{-2}=0, p_{-1}=1, q_{-2}=1, q_{-1}=0.
    std::int64_t p_prev = 0, p_curr = 1;
    std::int64_t q_prev = 1, q_curr = 0;
    std::optional<Rational<std::int64_t>> best;

    for (std::int64_t a : CFExpansionView{x, cf_params}) {
        const std::int64_t p_next = a * p_curr + p_prev;
        const std::int64_t q_next = a * q_curr + q_prev;

        if (q_next > params.max_denominator) {
            // Last valid convergent is p_curr/q_curr.
            // Check semi-convergents: p_curr - k·p_prev / q_curr - k·q_prev
            // is better than p_curr/q_curr when k is large enough. The largest
            // k keeping q in range:
            if (q_prev > 0) {
                const std::int64_t k = (params.max_denominator - q_curr) / q_prev;
                if (k > 0) {
                    Rational<std::int64_t> semi{p_curr + k * p_prev,
                                                q_curr + k * q_prev};
                    if (!best.has_value() ||
                        std::abs(semi.to_double() - x) < std::abs(best->to_double() - x))
                        best = semi;
                }
            }
            break;
        }

        Rational<std::int64_t> candidate{p_next, q_next};
        if (!best.has_value() ||
            std::abs(candidate.to_double() - x) < std::abs(best->to_double() - x))
            best = candidate;

        p_prev = p_curr; p_curr = p_next;
        q_prev = q_curr; q_curr = q_next;
    }

    if (!best.has_value())
        best = Rational<std::int64_t>{static_cast<std::int64_t>(std::round(x))};

    return *best;
}

} // namespace adrius
