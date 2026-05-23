// Copyright (c) 2025 InGenifold Research LLC. MIT License.
#pragma once

#include <adrius/backend/eigen.hpp>

namespace adrius {

// The backend used when no explicit Backend template argument is provided.
// Override by specializing adrius::DefaultBackend in your own code before
// including any algorithm header, or by passing an explicit Backend type.
using DefaultBackend = EigenBackend;

} // namespace adrius
