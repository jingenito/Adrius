# Boost.Multiprecision Backend

This document explains how to use the Boost.Multiprecision backend for arbitrary-precision arithmetic in Adrius.

## Overview

The default Eigen backend uses 64-bit floating-point (`double`), which provides approximately 15-17 significant decimal digits of precision. This is sufficient for most applications, but breaks down for:

- **Ill-conditioned problems**: Matrices with very different column magnitudes
- **High-denominator approximations**: Finding approximations with q > 10^10 requires higher precision
- **Very tight error bounds**: Approximation quality below 1e-15 requires higher precision

The **Boost.Multiprecision backend** allows arbitrary precision computation, with precision configurable from 50 to 1000+ decimal digits.

## Building with Boost.Multiprecision

### Option 1: Use System Boost
If Boost 1.70+ is already installed:
```bash
cmake -B build -DADRIUS_ENABLE_BOOST_MULTIPRECISION=ON -DADRIUS_USE_SYSTEM_BOOST=ON
cmake --build build
```

### Option 2: Auto-fetch Boost (CMake FetchContent)
CMake will download and build Boost automatically:
```bash
cmake -B build -DADRIUS_ENABLE_BOOST_MULTIPRECISION=ON
cmake --build build
# Note: First build may take 2-5 minutes to build Boost
```

### Disable Again
To switch back to Eigen-only builds:
```bash
cmake -B build -DADRIUS_ENABLE_BOOST_MULTIPRECISION=OFF
```

## Usage

### Basic Usage: Default Precision (50 digits)

```cpp
#include <adrius/adrius.hpp>
#include <iostream>
#include <numbers>
#include <span>

int main() {
    // Use Boost backend with default 50 decimal digits
    const std::vector<adrius::BoostBackendDefault::scalar_type> alpha = {
        adrius::BoostBackendDefault::scalar_type(std::numbers::sqrt2),
        adrius::BoostBackendDefault::scalar_type(std::numbers::pi)
    };

    adrius::ILLLParams params;
    params.max_iterations  = 30;
    params.max_denominator = 1'000'000'000LL;  // 1 billion
    
    auto result = adrius::illl<adrius::BoostBackendDefault>(
        std::span<const adrius::BoostBackendDefault::scalar_type>{alpha},
        params
    );

    for (std::size_t k = 0; k < result.iterations; ++k) {
        std::cout << "q = " << result.relations[k][0] << '\n';
    }
    
    return 0;
}
```

### Custom Precision

Change precision by specifying the `Digits` template parameter:

```cpp
// 100 decimal digits (332 bits)
using HighPrecision = adrius::BoostBackend<100>;
auto result = adrius::illl<HighPrecision>(alpha, params);

// 200 decimal digits (664 bits)
using VeryHighPrecision = adrius::BoostBackend<200>;
auto result = adrius::illl<VeryHighPrecision>(alpha, params);
```

**Predefined aliases:**
- `BoostBackendDefault` — 50 digits (recommended for most uses)
- `BoostBackendHighPrecision` — 100 digits
- `BoostBackendVeryHighPrecision` — 200 digits

### Type Conversion

Converting between Boost multiprecision and built-in types:

```cpp
using Backend = adrius::BoostBackendDefault;

// String → Boost scalar
auto x = Backend::scalar_type("3.141592653589793238462643383279502884197");

// Boost scalar → double (may lose precision)
double val = adrius::BoostBackend<50>::to_double(x);

// int → Boost scalar
Backend::scalar_type y = Backend::from_integer(42);

// Boost scalar → integer (truncates)
auto q = Backend::to_integer(x);  // Returns 3
```

### Comparing Backends

```cpp
#include <adrius/adrius.hpp>
#include <chrono>
#include <iostream>

template <typename Backend>
void benchmark_illl(const char* name, auto& alpha) {
    auto start = std::chrono::high_resolution_clock::now();
    
    auto result = adrius::illl<Backend>(
        std::span<decltype(alpha[0])>{alpha},
        adrius::ILLLParams{.max_iterations = 20}
    );
    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << name << ": " << elapsed.count() << " ms, "
              << result.iterations << " iterations\n";
}

int main() {
    // Eigen backend (built-in)
    std::vector<double> alpha_eigen = {1.414213562373095, 3.141592653589793};
    benchmark_illl<adrius::EigenBackend>("Eigen", alpha_eigen);

    // Boost backend (50 digits)
    std::vector<adrius::BoostBackendDefault::scalar_type> alpha_boost = {
        adrius::BoostBackendDefault::scalar_type("1.414213562373095"),
        adrius::BoostBackendDefault::scalar_type("3.141592653589793")
    };
    benchmark_illl<adrius::BoostBackendDefault>("Boost (50 digits)", alpha_boost);
    
    return 0;
}
```

## Performance Characteristics

Precision vs. Performance (approximate on modern hardware):

| Precision | Digits | Bits | Speed (rel. to Eigen) | Use Case |
|-----------|--------|------|----------------------|----------|
| Double (Eigen) | 15-17 | 64 | 1.0× | Standard problems |
| 50 digits (Boost default) | 50 | 166 | 10-20× slower | Most ill-conditioned cases |
| 100 digits (Boost high) | 100 | 332 | 50-100× slower | Very high-denominator approximations |
| 200 digits (Boost very high) | 200 | 664 | 200-300× slower | Extreme precision requirements |

**Recommendation**: Start with the Eigen backend. If you need higher precision, switch to Boost with 50 digits. Increase only if you hit numerical precision limits.

## Error Handling

Both backends raise the same exception types:

```cpp
#include <adrius/adrius.hpp>

try {
    auto result = adrius::illl<adrius::BoostBackendDefault>(alpha, params);
} catch (const adrius::DomainError& e) {
    std::cerr << "Bad input: " << e.what() << '\n';
} catch (const adrius::ConvergenceError& e) {
    std::cerr << "Algorithm didn't converge: " << e.what() << '\n';
} catch (const adrius::PrecisionError& e) {
    std::cerr << "Precision insufficient: " << e.what() << '\n';
}
```

## Debugging with Different Precisions

To diagnose if a problem is due to floating-point precision:

```cpp
// Try same algorithm with increasing precision
std::vector<double> alpha = {sqrt(2), M_PI};

// Step 1: Eigen
auto eigen_result = adrius::illl<adrius::EigenBackend>(
    /* convert to double */, params);

// Step 2: Boost 50 digits
auto boost50_result = adrius::illl<adrius::BoostBackend<50>>(
    /* convert to boost */, params);

// Step 3: Boost 100 digits (if Step 2 differs significantly from Step 1)
auto boost100_result = adrius::illl<adrius::BoostBackend<100>>(
    /* convert to boost */, params);

// If results stabilize at higher precision, the problem is ill-conditioned
// relative to double precision (not an algorithm bug).
```

## Implementation Details

The Boost backend uses:
- **`cpp_dec_float<Digits>`** for floating-point scalars (decimal, arbitrary precision)
- **`cpp_int`** for integer scalars (unbounded precision)
- **Eigen::Matrix** for matrix operations (works with any scalar type)

All standard Boost multiprecision functions (`sqrt`, `pow`, `abs`, etc.) are available via the Backend interface and match the Eigen API.

## References

- **Boost.Multiprecision Documentation**: https://www.boost.org/doc/libs/latest/libs/multiprecision/doc/html/index.html
- **cpp_dec_float**: https://www.boost.org/doc/libs/latest/libs/multiprecision/doc/html/boost_multiprecision/tut/floats/cpp_dec_float.html
- **Building Boost**: https://www.boost.org/doc/libs/latest/more/getting_started/

## See Also

- **DEVELOPMENT.md** — Build instructions and troubleshooting
- **CLAUDE.md** — Backend abstraction rules and design
- **docs/design-decisions.md** — Why Adrius uses concept-based backends
