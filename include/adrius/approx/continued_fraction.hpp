// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/core/error.hpp>
#include <adrius/core/params.hpp>

#include <cmath>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <vector>

namespace adrius {

// ── Lazy continued-fraction range ───────────────────────────────────────────
//
// CFExpansionView is a C++20 input range that yields partial quotients
// [a_0; a_1, a_2, …] on demand, terminating when max_depth is reached or
// the fractional part falls below zero_threshold.
//
// This is preferred over cf_expansion() when the caller wants early
// termination, lazy composition with std::views, or memory-bounded iteration.
//
// Example:
//   for (auto a : adrius::cf_view(std::numbers::phi) | std::views::take(10))
//       std::cout << a << ' ';
class CFExpansionView : public std::ranges::view_interface<CFExpansionView> {
public:
    struct Sentinel {};

    class Iterator {
    public:
        using value_type        = std::int64_t;
        using difference_type   = std::ptrdiff_t;
        using iterator_concept  = std::input_iterator_tag;

        Iterator() noexcept : done_(true) {}

        explicit Iterator(double x, CFParams params) noexcept
            : x_(x)
            , remaining_(params.max_depth)
            , threshold_(params.zero_threshold)
            , done_(params.max_depth == 0)
        {
            if (!done_) step_first();
        }

        [[nodiscard]] value_type operator*() const noexcept { return a_; }

        Iterator& operator++() noexcept {
            const double frac = x_ - static_cast<double>(a_);
            if (remaining_ == 1 || frac < threshold_) {
                done_ = true;
            } else {
                x_ = 1.0 / frac;
                --remaining_;
                a_ = static_cast<value_type>(std::floor(x_));
            }
            return *this;
        }

        // post-increment required by input_iterator
        void operator++(int) noexcept { ++(*this); }

        [[nodiscard]] bool operator==(const Sentinel&) const noexcept { return done_; }

    private:
        double      x_{};
        value_type  a_{};
        std::size_t remaining_{};
        double      threshold_{};
        bool        done_{};

        void step_first() noexcept {
            a_ = static_cast<value_type>(std::floor(x_));
        }
    };

    explicit CFExpansionView(double x, CFParams params = {}) noexcept
        : x_(x), params_(params) {}

    [[nodiscard]] Iterator begin() const noexcept { return Iterator{x_, params_}; }
    [[nodiscard]] Sentinel end()   const noexcept { return {}; }

private:
    double   x_;
    CFParams params_;
};

static_assert(std::ranges::input_range<CFExpansionView>);

// Convenience factory — name mirrors std::views:: convention.
[[nodiscard]] inline CFExpansionView
cf_view(double x, CFParams params = {}) noexcept
{
    return CFExpansionView{x, params};
}

// ── Eager variant ─────────────────────────────────────────────────────────
//
// Returns all partial quotients as a vector. Prefer cf_view() when you only
// need a prefix or want to compose with range adaptors.
[[nodiscard]] inline std::vector<std::int64_t>
cf_expansion(double x, CFParams params = {})
{
    std::vector<std::int64_t> result;
    result.reserve(params.max_depth);
    for (std::int64_t a : CFExpansionView{x, params})
        result.push_back(a);
    return result;
}

// ── Convergents ───────────────────────────────────────────────────────────
//
// Struct holding one convergent p_k/q_k and the next convergent's denominator
// (needed by best_rational to check the Farey neighbour).
struct Convergent {
    std::int64_t p{};  // numerator
    std::int64_t q{};  // denominator (always > 0)
};

// Returns convergents [p_0/q_0, p_1/q_1, …] for the CF expansion of x.
// The recurrence is:
//   p_{-1} = 1, p_0 = a_0, p_k = a_k·p_{k-1} + p_{k-2}
//   q_{-1} = 0, q_0 = 1,   q_k = a_k·q_{k-1} + q_{k-2}
[[nodiscard]] inline std::vector<Convergent>
cf_convergents(double x, CFParams params = {})
{
    std::vector<Convergent> result;
    result.reserve(params.max_depth);

    // Standard convergent recurrence seed: p_{-2}=0, p_{-1}=1, q_{-2}=1, q_{-1}=0.
    // First step gives p_0 = a_0·1 + 0 = a_0, q_0 = a_0·0 + 1 = 1, as required.
    std::int64_t p_prev = 0, p_curr = 1;
    std::int64_t q_prev = 1, q_curr = 0;

    for (std::int64_t a : CFExpansionView{x, params}) {
        const std::int64_t p_next = a * p_curr + p_prev;
        const std::int64_t q_next = a * q_curr + q_prev;
        result.push_back({p_next, q_next});
        p_prev = p_curr; p_curr = p_next;
        q_prev = q_curr; q_curr = q_next;
    }

    return result;
}

} // namespace adrius
