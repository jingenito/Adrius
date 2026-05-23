// Copyright (c) 2025 InGenifold Research LLC. MIT License.

#include <adrius/approx/continued_fraction.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include <ranges>
#include <vector>

TEST(ContinuedFraction, IntegerInputIsOneQuotient) {
    auto qs = adrius::cf_expansion(5.0);
    ASSERT_EQ(qs.size(), 1u);
    EXPECT_EQ(qs[0], 5);
}

TEST(ContinuedFraction, HalfGivesZeroOne) {
    // 0.5 = [0; 2]
    auto qs = adrius::cf_expansion(0.5);
    ASSERT_EQ(qs.size(), 2u);
    EXPECT_EQ(qs[0], 0);
    EXPECT_EQ(qs[1], 2);
}

TEST(ContinuedFraction, SimpleRational) {
    // 1.5 is exactly representable in binary (= 3/2 = [1; 2]).
    // 7/3 cannot be used here: the floating-point rounding of 7.0/3.0 causes
    // the second quotient to floor to 2 instead of 3, yielding [2; 2, 1].
    auto qs = adrius::cf_expansion(1.5);
    ASSERT_GE(qs.size(), 2u);
    EXPECT_EQ(qs[0], 1);
    EXPECT_EQ(qs[1], 2);
}

TEST(ContinuedFraction, GoldenRatioIsAllOnes) {
    // φ = [1; 1, 1, 1, …] — slowest converging CF
    auto qs = adrius::cf_expansion(std::numbers::phi, {.max_depth = 20});
    EXPECT_EQ(qs[0], 1);
    for (std::size_t i = 1; i < qs.size(); ++i)
        EXPECT_EQ(qs[i], 1) << "Quotient " << i << " should be 1 for the golden ratio";
}

TEST(ContinuedFraction, MaxDepthRespected) {
    adrius::CFParams p{.max_depth = 5};
    auto qs = adrius::cf_expansion(std::numbers::pi, p);
    EXPECT_LE(qs.size(), 5u);
}

TEST(ContinuedFraction, LazyViewComposesWithRanges) {
    // Use cf_view() + std::views::take — should compute only 4 quotients.
    // std::vector's iterator-pair ctor requires LegacyInputIterator; the
    // counted_iterator produced by views::take does not qualify on MSVC.
    // Use ranges::copy + back_inserter instead.
    auto view = adrius::cf_view(std::numbers::e) | std::views::take(4);
    std::vector<std::int64_t> qs;
    qs.reserve(4);
    std::ranges::copy(view, std::back_inserter(qs));
    ASSERT_EQ(qs.size(), 4u);
    // e = [2; 1, 2, 1, 1, 4, 1, 1, 6, …]
    EXPECT_EQ(qs[0], 2);
    EXPECT_EQ(qs[1], 1);
    EXPECT_EQ(qs[2], 2);
    EXPECT_EQ(qs[3], 1);
}

TEST(ContinuedFraction, ConvergentsRecurrenceHolds) {
    // For every convergent p_k/q_k: p_k·q_{k-1} - q_k·p_{k-1} = ±1
    auto convs = adrius::cf_convergents(std::numbers::pi, {.max_depth = 10});
    ASSERT_GE(convs.size(), 2u);
    for (std::size_t k = 1; k < convs.size(); ++k) {
        const auto cross = convs[k].p * convs[k-1].q - convs[k].q * convs[k-1].p;
        EXPECT_TRUE(cross == 1 || cross == -1)
            << "Convergent identity failed at k=" << k << " (cross=" << cross << ")";
    }
}
