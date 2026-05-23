// Copyright (c) 2025 InGenifold Research LLC. MIT License.

// Tests for BoostBackend: concept compliance, interface parity with EigenBackend,
// high-precision constants verified against published values, and ILLL integration.

#include <adrius/backend/boost_multiprecision.hpp>
#include <adrius/approx/illl.hpp>

#include <gtest/gtest.h>
#include <cstddef>
#include <numbers>
#include <span>

// ── Compile-time concept checks ───────────────────────────────────────────────

static_assert(adrius::Backend<adrius::BoostBackend<50>>,
    "BoostBackend<50> must satisfy the Backend concept");
static_assert(adrius::Backend<adrius::BoostBackend<100>>,
    "BoostBackend<100> must satisfy the Backend concept");
static_assert(adrius::Backend<adrius::BoostBackend<200>>,
    "BoostBackend<200> must satisfy the Backend concept");

// ── Type aliases ──────────────────────────────────────────────────────────────

using B50  = adrius::BoostBackend<50>;
using B100 = adrius::BoostBackend<100>;
using S50  = B50::scalar_type;
using S100 = B100::scalar_type;
using I    = B50::integer_type;

namespace {

// atan(x) via Gregory–Leibniz series for |x| < 1.
// Machin's formula π/4 = 4·atan(1/5) − atan(1/239) is used below for π.
template <typename S>
S atan_series(S x, int terms)
{
    S result{0}, x_sq = x * x, x_pow = x;
    for (int k = 0; k < terms; ++k) {
        S term = x_pow / S(2 * k + 1);
        result += (k % 2 == 0) ? term : -term;
        x_pow *= x_sq;
    }
    return result;
}

// e via Taylor series: sum_{k=0}^{terms} 1/k!
template <typename S>
S e_series(int terms)
{
    S result{1}, fact{1};
    for (int k = 1; k <= terms; ++k) {
        fact *= S(k);
        result += S(1) / fact;
    }
    return result;
}

} // namespace

// ── Backend interface tests (mirrors ConceptsEigenBackend) ────────────────────

TEST(ConceptsBoostBackend, FactoryFunctions) {
    auto M  = B50::make_zero_matrix(3, 4);
    auto v  = B50::make_zero_vector(3);
    auto Id = B50::identity_int(3);

    EXPECT_EQ(B50::rows(M), 3u);
    EXPECT_EQ(B50::cols(M), 4u);
    EXPECT_EQ(B50::norm_sq(v), S50{0});

    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_EQ(B50::int_get(Id, i, j), (i == j ? I{1} : I{0}));
}

TEST(ConceptsBoostBackend, ElementAccess) {
    auto M = B50::make_zero_matrix(2, 2);
    B50::set(M, 0, 1, S50("3.14159265358979323846264338327950288419716939937510"));
    EXPECT_EQ(B50::get(M, 0, 1),
              S50("3.14159265358979323846264338327950288419716939937510"));
    EXPECT_EQ(B50::get(M, 1, 0), S50{0});
}

TEST(ConceptsBoostBackend, IntegerAccess) {
    auto U = B50::identity_int(2);
    B50::int_set(U, 0, 1, I{"999999999999999999999999"});
    EXPECT_EQ(B50::int_get(U, 0, 1), I{"999999999999999999999999"});
    EXPECT_EQ(B50::int_get(U, 1, 1), I{1});
}

TEST(ConceptsBoostBackend, ColumnAccess) {
    auto M = B50::make_zero_matrix(3, 2);
    B50::set(M, 0, 0, S50{1}); B50::set(M, 1, 1, S50{1});

    auto c0 = B50::col(M, 0);
    auto c1 = B50::col(M, 1);
    EXPECT_EQ(B50::inner_product(c0, c1), S50{0});
    EXPECT_EQ(B50::norm_sq(c0), S50{1});

    auto v = B50::make_zero_vector(3);
    B50::set_col(M, 1, v);
    EXPECT_EQ(B50::norm_sq(B50::col(M, 1)), S50{0});
}

TEST(ConceptsBoostBackend, VectorArithmetic) {
    auto M = B50::make_zero_matrix(2, 2);
    B50::set(M, 0, 0, S50{3}); B50::set(M, 1, 1, S50{4});

    auto u = B50::col(M, 0);  // [3, 0]
    auto v = B50::col(M, 1);  // [0, 4]

    EXPECT_EQ(B50::inner_product(u, v), S50{0});
    EXPECT_EQ(B50::norm_sq(u), S50{9});
    EXPECT_EQ(B50::norm_sq(v), S50{16});
    EXPECT_EQ(B50::norm_sq(B50::add(u, v)), S50{25});      // 3-4-5 triangle
    EXPECT_EQ(B50::norm_sq(B50::subtract(u, v)), S50{25});
    EXPECT_EQ(B50::norm_sq(B50::scale(u, S50{2})), S50{36});
}

// ── Precision: sqrt(2) ────────────────────────────────────────────────────────

// NIST/published: sqrt(2) to 50 decimal places
// Source: https://oeis.org/A002193
static constexpr const char* kSqrt2_50 =
    "1.41421356237309504880168872420969807856967187537694";

TEST(BoostPrecision, SqrtTwoSquaredIdentity) {
    // sqrt(2)^2 must equal 2 within the precision of the 50-digit backend.
    S50 s = boost::multiprecision::sqrt(S50{2});
    S50 err = boost::multiprecision::abs(s * s - S50{2});
    EXPECT_LT(err, S50("1e-48"));
}

TEST(BoostPrecision, SqrtTwoKnownDigits50) {
    S50 computed = boost::multiprecision::sqrt(S50{2});
    S50 known(kSqrt2_50);
    S50 err = boost::multiprecision::abs(computed - known);
    EXPECT_LT(err, S50("1e-48"));
}

TEST(BoostPrecision, SqrtTwoKnownDigits100) {
    // NIST/OEIS: sqrt(2) to 100 decimal places (A002193)
    S100 computed = boost::multiprecision::sqrt(S100{2});
    S100 known(
        "1.4142135623730950488016887242096980785696718753769"
        "4807317667973799073247846210703885038753432764157274");
    S100 err = boost::multiprecision::abs(computed - known);
    EXPECT_LT(err, S100("1e-98"));
}

TEST(BoostPrecision, HigherDigitCountIsMoreAccurate) {
    // Embed the 50-digit result in 100-digit arithmetic and verify it falls
    // short of the 100-digit result by roughly 50 digits.
    S100 s50(kSqrt2_50);
    S100 s100 = boost::multiprecision::sqrt(S100{2});
    S100 diff = boost::multiprecision::abs(s100 - s50);
    // The 50-digit approximation should differ from the 100-digit one by ~1e-50
    EXPECT_GT(diff, S100("1e-55"));
    EXPECT_LT(diff, S100("1e-45"));
}

// ── Known constants ───────────────────────────────────────────────────────────

TEST(BoostConstants, PiViaMachinFormula50Digits) {
    // Machin (1706): π/4 = 4·atan(1/5) − atan(1/239)
    // Converges: 45 terms of atan(1/5) and 15 of atan(1/239) give >50 digits.
    // Published value (OEIS A000796):
    static constexpr const char* kPi50 =
        "3.14159265358979323846264338327950288419716939937510";

    S50 pi = S50{4} * (S50{4} * atan_series(S50{1}/S50{5}, 45)
                      - atan_series(S50{1}/S50{239}, 15));
    S50 err = boost::multiprecision::abs(pi - S50(kPi50));
    EXPECT_LT(err, S50("1e-48"));
}

TEST(BoostConstants, PiViaMachinFormula100Digits) {
    // Same formula at 100 digits; needs ~90 and ~30 terms respectively.
    // Published value (OEIS A000796):
    S100 kPi100(
        "3.1415926535897932384626433832795028841971693993751"
        "05820974944592307816406286208998628034825342117067982");

    S100 pi = S100{4} * (S100{4} * atan_series(S100{1}/S100{5}, 90)
                        - atan_series(S100{1}/S100{239}, 30));
    S100 err = boost::multiprecision::abs(pi - kPi100);
    EXPECT_LT(err, S100("1e-98"));
}

TEST(BoostConstants, EViaTaylorSeries50Digits) {
    // e = sum_{k=0}^inf 1/k!   — 70 terms exceeds 50 digits.
    // Published value (OEIS A001113):
    static constexpr const char* kE50 =
        "2.71828182845904523536028747135266249775724709369995";

    S50 e = e_series<S50>(70);
    S50 err = boost::multiprecision::abs(e - S50(kE50));
    EXPECT_LT(err, S50("1e-48"));
}

TEST(BoostConstants, GoldenRatioQuadraticIdentity) {
    // φ = (1 + sqrt(5)) / 2 satisfies φ² − φ − 1 = 0 exactly.
    S50 phi = (S50{1} + boost::multiprecision::sqrt(S50{5})) / S50{2};
    S50 residual = boost::multiprecision::abs(phi * phi - phi - S50{1});
    EXPECT_LT(residual, S50("1e-48"));
}

TEST(BoostConstants, GoldenRatioKnownDigits) {
    // Published value (OEIS A001622):
    static constexpr const char* kPhi50 =
        "1.61803398874989484820458683436563811772030917980576";
    S50 phi = (S50{1} + boost::multiprecision::sqrt(S50{5})) / S50{2};
    S50 err = boost::multiprecision::abs(phi - S50(kPhi50));
    EXPECT_LT(err, S50("1e-48"));
}

TEST(BoostConstants, PythagoreanIdentityHighPrecision) {
    // sin²(x) + cos²(x) = 1 for x = π/6 (using exact sin=1/2, cos=sqrt(3)/2).
    S50 sin_pi6 = S50{1} / S50{2};
    S50 cos_pi6 = boost::multiprecision::sqrt(S50{3}) / S50{2};
    S50 residual = boost::multiprecision::abs(sin_pi6*sin_pi6 + cos_pi6*cos_pi6 - S50{1});
    EXPECT_LT(residual, S50("1e-48"));
}

// ── Integer arithmetic ────────────────────────────────────────────────────────

TEST(BoostInteger, FiftyFactorialExact) {
    // 50! verified against published value (65 digits, 12 trailing zeros).
    // Source: NIST Handbook of Mathematical Functions, Table 5.1
    I result{1};
    for (int k = 2; k <= 50; ++k) result *= I{k};
    EXPECT_EQ(result,
        I{"30414093201713378043612608166979581188299763898377856000000000000"});
}

TEST(BoostInteger, TrailingZerosOf50Factorial) {
    // Legendre's formula: trailing zeros of n! = sum floor(n/5^k)
    // 50: 10 + 2 = 12 trailing zeros.
    I result{1};
    for (int k = 2; k <= 50; ++k) result *= I{k};
    I trailing{0};
    I tmp = result;
    while (tmp % I{10} == I{0}) { ++trailing; tmp /= I{10}; }
    EXPECT_EQ(trailing, I{12});
}

TEST(BoostInteger, LargeExponentExact) {
    // 2^100 = 1267650600228229401496703205376  (31 digits, verified)
    I result{1};
    for (int k = 0; k < 100; ++k) result *= I{2};
    EXPECT_EQ(result, I{"1267650600228229401496703205376"});
}

// ── ILLL integration with BoostBackend ───────────────────────────────────────

TEST(BoostILLL, Sqrt2ConvergentsMatchKnownCF) {
    // CF convergents of sqrt(2): 1/1, 3/2, 7/5, 17/12, 41/29, 99/70, 239/169
    // ILLL in 1D with denominator cap 300 should find (q,p) matching these.
    using B = B50;
    const std::vector<S50> alpha = {boost::multiprecision::sqrt(S50{2})};
    adrius::ILLLParams params;
    params.max_iterations  = 20;
    params.max_denominator = 300;
    params.quality_tol     = 0;

    auto result = adrius::illl<B>(std::span<const S50>{alpha}, params);
    ASSERT_GT(result.iterations, 0u);

    // Every relation must satisfy the Hurwitz bound: q·|q·√2 − p| ≤ 1/√5 < 0.448
    for (std::size_t k = 0; k < result.iterations; ++k) {
        S50 q   = S50(result.relations[k][0]);
        S50 err = result.quality[k];
        EXPECT_GT(q, S50{0});
        EXPECT_LT(q * err, S50{"0.45"})
            << "Hurwitz bound violated at step " << k;
    }

    // At least one of the known convergents must appear.
    // Convergent pairs (q, p): (1,1), (2,3), (5,7), (12,17), (29,41), (70,99), (169,239)
    using I50 = B50::integer_type;
    const std::vector<std::pair<I50,I50>> known_convergents = {
        {1,1},{2,3},{5,7},{12,17},{29,41},{70,99},{169,239}
    };
    bool found = false;
    for (const auto& rel : result.relations) {
        I50 q = rel[0], p = rel[1];
        for (const auto& [kq, kp] : known_convergents)
            if (q == kq && p == kp) { found = true; break; }
        if (found) break;
    }
    EXPECT_TRUE(found) << "No known CF convergent of sqrt(2) found";
}

TEST(BoostILLL, DirichletBoundHoldsForPiAndE) {
    // 2D approximation: alpha = (pi, e) — a genuinely hard simultaneous case.
    // Bosma-Smeets guarantees for n=2: q^{1/2}·max|q·αᵢ − pᵢ| ≤ C(2).
    using B = B50;
    adrius::ILLLParams params;
    params.max_iterations  = 15;
    params.max_denominator = 0;
    params.quality_tol     = 0;

    // pi via Machin to 50 digits, e via Taylor
    S50 pi50 = S50{4} * (S50{4} * atan_series(S50{1}/S50{5}, 45)
                        - atan_series(S50{1}/S50{239}, 15));
    S50 e50  = e_series<S50>(70);

    const std::vector<S50> alpha = {pi50, e50};
    auto result = adrius::illl<B>(std::span<const S50>{alpha}, params);
    ASSERT_GT(result.iterations, 0u);

    for (std::size_t k = 0; k < result.iterations; ++k) {
        S50 q = S50(result.relations[k][0]);
        EXPECT_GT(q, S50{0}) << "Zero denominator at step " << k;
        // Generous bound: q^{1/2} · max_err < 2.0
        S50 dirichlet = boost::multiprecision::sqrt(q) * result.quality[k];
        EXPECT_LT(dirichlet, S50{2})
            << "Dirichlet coefficient out of range at step " << k;
    }
}

TEST(BoostILLL, HighPrecisionFindsLargerDenominators) {
    // With 50-digit floats we can safely seek denominators far beyond what
    // double precision (15-17 digits) can represent without rounding errors.
    using B = B50;
    S50 sqrt2 = boost::multiprecision::sqrt(S50{2});
    const std::vector<S50> alpha = {sqrt2};

    adrius::ILLLParams params;
    params.max_iterations  = 40;
    params.max_denominator = 1000000;  // 10^6 — well beyond double's safe range for this
    params.quality_tol     = 0;

    auto result = adrius::illl<B>(std::span<const S50>{alpha}, params);
    ASSERT_GT(result.iterations, 0u);

    // The largest denominator found should substantially exceed what
    // a 4-iteration double-precision ILLL would reach (typically < 1000).
    I max_q = result.relations[0][0];
    for (const auto& rel : result.relations)
        if (rel[0] > max_q) max_q = rel[0];

    EXPECT_GT(max_q, I{1000}) << "Expected large denominators with 50-digit precision";
}
