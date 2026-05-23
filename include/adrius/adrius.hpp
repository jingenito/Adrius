// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

// Umbrella header — include this for the full library with the default
// (Eigen3) backend.
//
// Selective inclusion: if you are using a custom backend, include only the
// algorithm headers you need (e.g. <adrius/linalg/lll.hpp>) plus your own
// backend header. Do NOT include this file in that case.

#include <adrius/version.hpp>

// Core (no external dependencies)
#include <adrius/core/concepts.hpp>
#include <adrius/core/error.hpp>
#include <adrius/core/params.hpp>
#include <adrius/core/traits.hpp>

// Backend (pulls in Eigen3; optionally Boost.Multiprecision)
#include <adrius/backend/eigen.hpp>
#include <adrius/backend/default_backend.hpp>
#ifdef ADRIUS_HAS_BOOST_MULTIPRECISION
#include <adrius/backend/boost_multiprecision.hpp>
#endif

// Utilities
#include <adrius/util/rational_type.hpp>

// Linear algebra
#include <adrius/linalg/gram_schmidt.hpp>
#include <adrius/linalg/lll.hpp>

// Approximation algorithms
#include <adrius/approx/continued_fraction.hpp>
#include <adrius/approx/rational.hpp>
#include <adrius/approx/simultaneous.hpp>
#include <adrius/approx/illl.hpp>
