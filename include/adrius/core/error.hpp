// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <stdexcept>

namespace adrius {

// Root of the adrius exception hierarchy.
struct Error : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// A mathematical precondition was violated — linearly dependent basis,
// zero denominator, parameter outside its valid range, etc.
struct DomainError : Error {
    using Error::Error;
};

// An iterative algorithm exhausted its iteration budget without converging.
struct ConvergenceError : Error {
    using Error::Error;
};

// Floating-point precision is insufficient for the magnitude or conditioning
// of the given input.
struct PrecisionError : Error {
    using Error::Error;
};

} // namespace adrius
