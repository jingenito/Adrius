// Copyright (c) 2025 InGenifold Research LLC. MIT License.

// Verify the Backend concept is satisfied by EigenBackend at compile time,
// and exercise every required static member at runtime to catch linker gaps.

#include <adrius/core/concepts.hpp>
#include <adrius/backend/eigen.hpp>

#include <gtest/gtest.h>
#include <cstdint>

// Compile-time contract check — caught during compilation, not at test runtime.
static_assert(adrius::Backend<adrius::EigenBackend>,
    "EigenBackend must satisfy the Backend concept");

using B = adrius::EigenBackend;

TEST(ConceptsEigenBackend, TypeAliasesExist) {
    static_assert(std::floating_point<B::scalar_type>);
    static_assert(std::signed_integral<B::integer_type>);
    SUCCEED();
}

TEST(ConceptsEigenBackend, FactoryFunctions) {
    auto M  = B::make_zero_matrix(3, 4);
    auto v  = B::make_zero_vector(3);
    auto Id = B::identity_int(3);

    EXPECT_EQ(B::rows(M), 3u);
    EXPECT_EQ(B::cols(M), 4u);
    EXPECT_DOUBLE_EQ(B::norm_sq(v), 0.0);

    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            EXPECT_EQ(B::int_get(Id, i, j), (i == j ? std::int64_t{1} : std::int64_t{0}));
}

TEST(ConceptsEigenBackend, ElementAccess) {
    auto M = B::make_zero_matrix(2, 2);
    B::set(M, 0, 1, 3.14);
    EXPECT_DOUBLE_EQ(B::get(M, 0, 1), 3.14);
    EXPECT_DOUBLE_EQ(B::get(M, 1, 0), 0.0);
}

TEST(ConceptsEigenBackend, IntElementAccess) {
    auto U = B::identity_int(2);
    B::int_set(U, 0, 1, std::int64_t{7});
    EXPECT_EQ(B::int_get(U, 0, 1), std::int64_t{7});
    EXPECT_EQ(B::int_get(U, 1, 1), std::int64_t{1});
}

TEST(ConceptsEigenBackend, VectorArithmetic) {
    // Build u = [1, 0, 0], v = [0, 1, 0] via column access on a matrix
    auto basis = B::make_zero_matrix(3, 2);
    B::set(basis, 0, 0, 1.0);
    B::set(basis, 1, 1, 1.0);

    auto col0 = B::col(basis, 0);
    auto col1 = B::col(basis, 1);

    EXPECT_DOUBLE_EQ(B::inner_product(col0, col1), 0.0);
    EXPECT_DOUBLE_EQ(B::norm_sq(col0), 1.0);

    auto sum  = B::add(col0, col1);
    auto diff = B::subtract(col0, col1);
    auto sc   = B::scale(col0, 5.0);

    EXPECT_DOUBLE_EQ(B::norm_sq(sum),  2.0);
    EXPECT_DOUBLE_EQ(B::norm_sq(diff), 2.0);
    EXPECT_DOUBLE_EQ(B::norm_sq(sc),   25.0);
}

TEST(ConceptsEigenBackend, MapVectorSpan) {
    const double data[3] = {1.0, 2.0, 3.0};
    auto view = B::map_vector(std::span<const double>{data, 3});
    EXPECT_DOUBLE_EQ(B::norm_sq(view), 14.0);
}
