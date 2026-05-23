// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/core/concepts.hpp>

namespace adrius {

// Type aliases extracted from a Backend.
// Use these in algorithm signatures instead of typename B::xxx to reduce
// noise in deeply-nested template expressions.

template <Backend B> using scalar_of     = typename B::scalar_type;
template <Backend B> using integer_of    = typename B::integer_type;
template <Backend B> using matrix_of     = typename B::matrix_type;
template <Backend B> using vector_of     = typename B::vector_type;
template <Backend B> using int_matrix_of = typename B::int_matrix_type;

} // namespace adrius
