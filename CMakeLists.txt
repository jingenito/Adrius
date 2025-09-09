cmake_minimum_required(VERSION 3.20)

message(STATUS "=== STARTING CMAKE CONFIGURATION ===")

# Set vcpkg toolchain before project declaration
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
        CACHE STRING "")
endif()

project(adrius)

message(STATUS "=== PROJECT DECLARED ===")

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Check for C++20 support
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "GCC 10.0 or later is required for C++20 support")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Clang 10.0 or later is required for C++20 support")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.29")
    message(FATAL_ERROR "MSVC 19.29 or later is required for C++20 support")
endif()

message(STATUS "=== ABOUT TO FIND EIGEN3 ===")
message(STATUS "CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}")
message(STATUS "VCPKG_INSTALLED_DIR: ${VCPKG_INSTALLED_DIR}")
message(STATUS "VCPKG_TARGET_TRIPLET: ${VCPKG_TARGET_TRIPLET}")

# Find Eigen3 - use CONFIG mode for vcpkg
find_package(Eigen3 CONFIG REQUIRED)

message(STATUS "=== AFTER FIND_PACKAGE ===")

# Check if target exists
if(TARGET Eigen3::Eigen)
    message(STATUS "SUCCESS: Eigen3::Eigen target found!")
    get_target_property(EIGEN_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "Eigen3 include directories: ${EIGEN_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "FAILED: Eigen3::Eigen target not found!")
endif()

# Simple executable for testing
add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE Eigen3::Eigen)

message(STATUS "=== CMAKE CONFIGURATION COMPLETE ===")