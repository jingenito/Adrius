# Copyright (c) 2025 InGenifold Research LLC. MIT License.
cmake_minimum_required(VERSION 3.21)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# pybind11 — prefer a system/pip install, fall back to FetchContent.
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found via find_package; fetching from GitHub")
    include(FetchContent)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.13.6
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(pybind11)
endif()

# Build the C extension as _adrius.pyd / _adrius.so inside the adrius/ package
# directory so `from ._adrius import ...` works without any PYTHONPATH gymnastics.
pybind11_add_module(_adrius adrius_py.cpp)
target_link_libraries(_adrius PRIVATE adrius::adrius)
target_compile_options(_adrius PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)

# Force the .pyd / .so into build/python/adrius/ (works for all generators,
# including MSVC multi-config which would otherwise add a Debug/ subdirectory).
set(_pkg_dir "${CMAKE_CURRENT_BINARY_DIR}/adrius")
foreach(_prop
        LIBRARY_OUTPUT_DIRECTORY         RUNTIME_OUTPUT_DIRECTORY
        LIBRARY_OUTPUT_DIRECTORY_DEBUG   RUNTIME_OUTPUT_DIRECTORY_DEBUG
        LIBRARY_OUTPUT_DIRECTORY_RELEASE RUNTIME_OUTPUT_DIRECTORY_RELEASE
        LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO
        RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO)
    set_target_properties(_adrius PROPERTIES ${_prop} "${_pkg_dir}")
endforeach()

# Copy the pure-Python package source into the build tree alongside the extension.
file(COPY adrius/ DESTINATION "${_pkg_dir}")
