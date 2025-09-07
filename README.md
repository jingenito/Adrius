# Adrius

A modern C++ library for Diophantine approximation with pluggable linear algebra backends.

## Overview

Adrius provides efficient algorithms for Diophantine approximation problems, including rational approximations to real numbers, simultaneous approximations, and related number-theoretic computations. Unlike other libraries in this domain that require GNU licensing, Adrius is released under the permissive MIT license.

## Key Features

- **MIT Licensed**: Use freely in both open-source and commercial projects
- **Pluggable Backends**: Abstract linear algebra interface supports multiple libraries
- **High Performance**: Optimized algorithms with efficient matrix operations
- **Modern C++**: Clean, type-safe API with C++17/20 features
- **Header-Only**: Easy integration with CMake FetchContent or package managers

## Algorithms Supported

- Continued fraction expansions
- Best rational approximations
- Simultaneous Diophantine approximation
- Lattice reduction algorithms
- Padé approximation
- Convergents and semiconvergents

## Quick Start

### Installation

#### Using CMake FetchContent
```cmake
include(FetchContent)
FetchContent_Declare(
    adrius
    GIT_REPOSITORY https://github.com/jingenito/adrius
    GIT_TAG main
)
FetchContent_MakeAvailable(adrius)

target_link_libraries(your_target adrius::adrius)
```

#### Manual Installation
```bash
git clone https://github.com/jingenito/adrius.git
cd adrius
mkdir build && cd build
cmake ..
make install
```

### Basic Usage

```cpp
#include <adrius/adrius.hpp>

int main() {
    // Find best rational approximation to π with denominator ≤ 1000
    auto approx = adrius::best_approximation(M_PI, 1000);
    std::cout << "π ≈ " << approx.numerator << "/" << approx.denominator << std::endl;
    
    // Continued fraction expansion
    auto cf = adrius::continued_fraction(std::sqrt(2), 10);
    
    // Simultaneous approximation
    std::vector<double> targets = {M_PI, M_E, std::sqrt(2)};
    auto simultaneous = adrius::simultaneous_approximation(targets, 1000);
    
    return 0;
}
```

## Linear Algebra Backends

Adrius abstracts linear algebra operations through a pluggable interface. The default implementation uses Eigen, but you can easily switch backends:

### Default (Eigen)
```cpp
#include <adrius/backends/eigen_backend.hpp>
using Backend = adrius::EigenBackend<double>;
```

### Custom Backend
```cpp
#include <adrius/backends/custom_backend.hpp>
using Backend = adrius::CustomBackend<double>;
```

### Supported Backends
- **Eigen** (default) - Header-only, high performance
- **BLAS/LAPACK** - For maximum performance on large problems
- **Custom** - Implement your own for specialized hardware

## Configuration Options

```cpp
// Precision control
adrius::Config config;
config.precision = 1e-12;
config.max_iterations = 1000;

// Algorithm selection
config.algorithm = adrius::Algorithm::LLL;  // or PSLQ, Ferguson-Forcade
```

## Examples

See the `examples/` directory for complete working examples:
- `basic_approximation.cpp` - Simple rational approximations
- `continued_fractions.cpp` - Working with continued fractions
- `simultaneous.cpp` - Multi-dimensional approximation
- `custom_backend.cpp` - Using different linear algebra backends

## Dependencies

### Required
- C++17 compatible compiler
- CMake 3.12+

### Default Backend (Eigen)
- Eigen3 (automatically fetched if not found)

### Optional
- BLAS/LAPACK (for alternative backends)
- Boost (for extended precision types)

## Building from Source

```bash
git clone https://github.com/jingenito/adrius.git
cd adrius

# Configure
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build

# Run tests
cmake --build build --target test

# Install
cmake --install build
```

### CMake Options
- `ADRIUS_BUILD_TESTS=ON` - Build test suite
- `ADRIUS_BUILD_EXAMPLES=ON` - Build example programs  
- `ADRIUS_USE_SYSTEM_EIGEN=OFF` - Use system Eigen instead of fetching
- `ADRIUS_ENABLE_BENCHMARKS=OFF` - Build performance benchmarks

## Performance

Adrius is designed for high performance:
- Template-based design eliminates virtual function overhead
- Efficient algorithms with optimal complexity
- SIMD optimizations where applicable
- Memory pool allocation for frequent operations

Benchmark results on modern hardware show competitive performance with specialized libraries while maintaining the flexibility of pluggable backends.

## License Advantage

Most Diophantine approximation libraries require GNU GPL licensing, which restricts use in proprietary software. Adrius uses only MIT-compatible dependencies:

- **Eigen**: Mozilla Public License 2.0
- **Standard Library**: Implementation specific (typically permissive)

This makes Adrius suitable for:
- Commercial software development
- Embedded systems
- Academic research with publication requirements
- Integration into larger proprietary codebases

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup
```bash
git clone https://github.com/jingenito/adrius.git
cd adrius
cmake -B build -DADRIUS_BUILD_TESTS=ON -DADRIUS_BUILD_EXAMPLES=ON
cmake --build build
```

## Documentation

- [API Reference](docs/api.md)
- [Algorithm Details](docs/algorithms.md)
- [Backend Development](docs/backends.md)
- [Performance Guide](docs/performance.md)

## Citation

If you use Adrius in academic work, please cite:

```bibtex
@software{adrius,
  title={Adrius: A C++ Library for Diophantine Approximation},
  author={Your Name},
  url={https://github.com/jingenito/adrius},
  year={2025}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

Copyright (c) 2025 InGenifold Research LLC

## Acknowledgments

- The Eigen team for providing an excellent linear algebra library
- Contributors to the field of Diophantine approximation
- The C++ community for modern language features that make this library possible