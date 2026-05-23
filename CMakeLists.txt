cmake_minimum_required(VERSION 3.21)
project(adrius VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MSVC requires /EHsc for standard C++ exception handling. The Visual Studio
# generator injects this automatically; Ninja does not, so we set it explicitly.
if(MSVC)
    add_compile_options(/EHsc)
endif()

# Options
option(ADRIUS_BUILD_TESTS     "Build test suite"                       ON)
option(ADRIUS_BUILD_EXAMPLES  "Build example programs"                 ON)
option(ADRIUS_USE_SYSTEM_EIGEN "Use an installed Eigen3 (skip fetch)"  OFF)
# Guard install targets so that FetchContent consumers don't pollute their
# install prefix. Defaults ON only when this is the top-level project.
option(ADRIUS_ENABLE_INSTALL  "Generate install targets"
       ${PROJECT_IS_TOP_LEVEL})

include(FetchContent)

# --- Eigen3 ---
if(ADRIUS_USE_SYSTEM_EIGEN)
    find_package(Eigen3 3.4 REQUIRED NO_MODULE)
else()
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
        GIT_SHALLOW    TRUE
    )
    # Prevent Eigen from registering its own ~900 unit tests with CTest.
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen)
    # Restore so our own tests can re-enable testing below.
    set(BUILD_TESTING ON CACHE BOOL "" FORCE)
endif()

# --- Header-only interface target ---
add_library(adrius INTERFACE)
target_include_directories(adrius INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(adrius INTERFACE Eigen3::Eigen)
target_compile_features(adrius INTERFACE cxx_std_20)
add_library(adrius::adrius ALIAS adrius)

# --- Tests ---
if(ADRIUS_BUILD_TESTS)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.14.0
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    add_subdirectory(tests)
endif()

# --- Examples ---
if(ADRIUS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# --- Install / find_package support ---
if(ADRIUS_ENABLE_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(TARGETS adrius
        EXPORT AdriusTargets)

    install(DIRECTORY include/adrius
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(EXPORT AdriusTargets
        FILE        AdriusTargets.cmake
        NAMESPACE   adrius::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/adrius)

    configure_package_config_file(
        cmake/AdriusConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/AdriusConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/adrius)

    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/AdriusConfigVersion.cmake
        VERSION     ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)

    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/AdriusConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/AdriusConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/adrius)
endif()