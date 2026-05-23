// Copyright (c) 2025 InGenifold Research LLC. MIT License.

#include <adrius/util/rational_type.hpp>
#include <adrius/approx/rational.hpp>

#include <gtest/gtest.h>
#include <cmath>
#include <numbers>

using R = adrius::Rational<std::int64_t>;

// ── Rational<> arithmetic ────────────────────────────────────────────────────

TEST(Rational, DefaultIsZero) {
    R r;
    EXPECT_EQ(r.numerator(),   0);
    EXPECT_EQ(r.denominator(), 1);
    EXPECT_TRUE(r.is_zero());
}

TEST(Rational, ReducesToLowestTerms) {
    R r{6, 4};
    EXPECT_EQ(r.numerator(),   3);
    EXPECT_EQ(r.denominator(), 2);
}

TEST(Rational, NegativeDenominatorNormalized) {
    R r{3, -4};
    EXPECT_EQ(r.numerator(),   -3);
    EXPECT_EQ(r.denominator(),  4);
}

TEST(Rational, ZeroDenominatorThrows) {
    EXPECT_THROW((R{1, 0}), adrius::DomainError);
}

TEST(Rational, Addition) {
    EXPECT_EQ((R{1,2} + R{1,3}), R(5,6));
    EXPECT_EQ((R{1,2} + R{-1,2}), R(0));
}

TEST(Rational, Subtraction) {
    EXPECT_EQ((R{3,4} - R{1,4}), R(1,2));
}

TEST(Rational, Multiplication) {
    EXPECT_EQ((R{2,3} * R{3,4}), R(1,2));
}

TEST(Rational, Division) {
    EXPECT_EQ((R{1,2} / R{1,4}), R(2));
}

TEST(Rational, Negation) {
    EXPECT_EQ(-R(3,5), R(-3,5));
    EXPECT_EQ(-R(-3,5), R(3,5));
}

TEST(Rational, ComparisonOrdering) {
    EXPECT_LT(R(1,3), R(1,2));
    EXPECT_GT(R(2,3), R(1,2));
    EXPECT_EQ(R(2,4), R(1,2));
    EXPECT_LE(R(1,2), R(1,2));
}

TEST(Rational, ToDouble) {
    EXPECT_DOUBLE_EQ(R(1,4).to_double(), 0.25);
    EXPECT_DOUBLE_EQ(R(1,3).to_double(), 1.0/3.0);
}

TEST(Rational, IsInteger) {
    EXPECT_TRUE(R(6,3).is_integer());   // reduces to 2/1
    EXPECT_FALSE(R(5,3).is_integer());
}

// ── best_rational ─────────────────────────────────────────────────────────

TEST(BestRational, ExactInteger) {
    auto r = adrius::best_rational(3.0);
    EXPECT_EQ(r.numerator(),   3);
    EXPECT_EQ(r.denominator(), 1);
}

TEST(BestRational, ExactHalf) {
    auto r = adrius::best_rational(0.5);
    EXPECT_EQ(r.numerator(),   1);
    EXPECT_EQ(r.denominator(), 2);
}

TEST(BestRational, PiApproximation) {
    // 355/113 is the famous best rational approximation to π with q ≤ 1000.
    auto r = adrius::best_rational(std::numbers::pi, {.max_denominator = 1000});
    EXPECT_EQ(r.denominator(), 113);
    EXPECT_EQ(r.numerator(),   355);
}

TEST(BestRational, DenominatorRespected) {
    auto r = adrius::best_rational(std::numbers::pi, {.max_denominator = 10});
    EXPECT_LE(r.denominator(), 10);
}

TEST(BestRational, ZeroDenominatorThrows) {
    EXPECT_THROW((void)adrius::best_rational(1.0, {.max_denominator = 0}),
                 adrius::DomainError);
}
