# Adrius

A C++20 header-only library for Diophantine approximation.
MIT license — InGenifold Research LLC.

## Algorithms

| Header | What it provides |
|---|---|
| `adrius/linalg/gram_schmidt.hpp` | Classical Gram-Schmidt orthogonalization |
| `adrius/linalg/lll.hpp` | LLL lattice reduction (Lovász incremental updates) |
| `adrius/approx/continued_fraction.hpp` | Lazy CF expansion range + convergents |
| `adrius/approx/rational.hpp` | Best rational approximation via convergents |
| `adrius/approx/simultaneous.hpp` | Simultaneous Diophantine approximation (LLL-based) |
| `adrius/approx/illl.hpp` | ILLL — Iterated LLL (Bosma & Smeets 2010) |
| `adrius/util/rational_type.hpp` | `Rational<Int>` — exact reduced-form arithmetic |

Include `<adrius/adrius.hpp>` to pull in everything with the default Eigen3 backend.

## Quick start

### Continued fraction expansion

```cpp
#include <adrius/adrius.hpp>
#include <iostream>
#include <numbers>
#include <ranges>

int main() {
    // Lazy range — computes quotients on demand.
    // Compose freely with std::views:
    for (std::int64_t a : adrius::cf_view(std::numbers::pi) | std::views::take(8))
        std::cout << a << ' ';   // 3 7 15 1 292 1 1 1
    std::cout << '\n';

    // Eager version when you want a vector:
    auto qs = adrius::cf_expansion(std::numbers::phi, {.max_depth = 10});

    // Convergents p_k / q_k:
    auto convs = adrius::cf_convergents(std::numbers::pi, {.max_depth = 6});
    for (auto [p, q] : convs)
        std::cout << p << '/' << q << '\n';  // 3/1  22/7  333/106  ...
}
```

### Best rational approximation

```cpp
#include <adrius/adrius.hpp>
#include <iostream>
#include <numbers>

int main() {
    // 355/113 is the best approximation to π with denominator ≤ 1000.
    auto r = adrius::best_rational(std::numbers::pi, {.max_denominator = 1000});
    std::cout << r << '\n';  // 355/113
}
```

### Simultaneous Diophantine approximation

```cpp
#include <adrius/adrius.hpp>
#include <iostream>
#include <numbers>
#include <span>

int main() {
    // Find q, p₁, p₂ such that |q·αᵢ − pᵢ| is small for both αs.
    const std::vector<double> alpha = {std::numbers::sqrt2, std::numbers::pi};

    auto result = adrius::simultaneous_approx<adrius::EigenBackend>(
        std::span<const double>{alpha}, /*scale=*/1e6);

    std::cout << "q  = " << result.denominator << '\n';
    std::cout << "p  = " << result.numerators[0] << ", " << result.numerators[1] << '\n';
    std::cout << "err = " << result.quality << '\n';
}
```

### LLL lattice reduction

```cpp
#include <adrius/adrius.hpp>

int main() {
    // Build a 3×3 basis (columns are basis vectors).
    adrius::EigenBackend::matrix_type basis(3, 3);
    basis << 1, -1,  3,
             1,  0,  5,
             1,  2,  6;

    // Two-stage API: preprocess validates + computes initial GSO,
    // lll_reduce iterates with Lovász incremental updates.
    auto prepared = adrius::preprocess_lll<adrius::EigenBackend>(basis);
    auto result   = adrius::lll_reduce<adrius::EigenBackend>(std::move(prepared));

    // result.reduced_basis  — the LLL-reduced columns
    // result.transform      — integer matrix U: reduced = original · U
    // result.swap_count     — number of Lovász swaps performed
}
```

### Arbitrary-Precision Backend (Boost.Multiprecision)

For ill-conditioned problems or high-precision approximations, use the Boost.Multiprecision backend:

```cpp
#include <adrius/adrius.hpp>

int main() {
    // 50 decimal digits (default)
    using BoostPrecision = adrius::BoostBackendDefault;
    
    std::vector<BoostPrecision::scalar_type> alpha = {
        BoostPrecision::scalar_type("1.414213562373095"),
        BoostPrecision::scalar_type("3.141592653589793")
    };
    
    auto result = adrius::illl<BoostPrecision>(
        std::span<const BoostPrecision::scalar_type>{alpha},
        adrius::ILLLParams{.max_denominator = 1'000'000'000LL}
    );
}
```

Build with: `cmake -B build -DADRIUS_ENABLE_BOOST_MULTIPRECISION=ON`

See `docs/boost-multiprecision-backend.md` for detailed usage, precision options, and performance considerations.

### Custom Backend

Algorithms are written against the `adrius::Backend` concept — no Eigen dependency in algorithm headers.

```cpp
// my_backend.hpp
#include <adrius/core/concepts.hpp>

struct MyBackend {
    using scalar_type     = float;
    using integer_type    = std::int32_t;
    using matrix_type     = /* your type */;
    using vector_type     = /* your type */;
    using int_matrix_type = /* your type */;

    static std::size_t rows(const matrix_type&);
    static std::size_t cols(const matrix_type&);
    // ... (see docs/design-decisions.md §1 for the full required interface)
};

static_assert(adrius::Backend<MyBackend>);

// Then use it:
#include <adrius/linalg/lll.hpp>
auto result = adrius::lll_reduce<MyBackend>(basis);
```

## Building

Requires CMake 3.21+ and a C++20 compiler. Eigen3 and GoogleTest are fetched automatically.

```bash
# Configure (Ninja + compile_commands.json for IntelliSense)
cmake -B build --preset windows-msvc-debug   # Windows
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON     # Linux / macOS

# Build
cmake --build build

# Test
ctest --test-dir build --output-on-failure

# Install (supports find_package)
cmake --install build
```

### CMake options

| Option | Default | Description |
|---|---|---|
| `ADRIUS_BUILD_TESTS` | `ON` | Build GoogleTest suite |
| `ADRIUS_BUILD_EXAMPLES` | `ON` | Build example programs |
| `ADRIUS_USE_SYSTEM_EIGEN` | `OFF` | Use an installed Eigen3 instead of fetching |
| `ADRIUS_ENABLE_BOOST_MULTIPRECISION` | `OFF` | Enable Boost.Multiprecision arbitrary-precision backend |
| `ADRIUS_USE_SYSTEM_BOOST` | `OFF` | Use an installed Boost instead of fetching (requires `ADRIUS_ENABLE_BOOST_MULTIPRECISION=ON`) |
| `ADRIUS_ENABLE_INSTALL` | `ON` (top-level) | Generate install/find_package targets |

### Consuming via FetchContent

```cmake
include(FetchContent)
FetchContent_Declare(adrius
    GIT_REPOSITORY https://github.com/jingenito/Adrius.git
    GIT_TAG        main
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(adrius)

target_link_libraries(your_target PRIVATE adrius::adrius)
```

### Consuming via find_package (after install)

```cmake
find_package(adrius 0.1 REQUIRED)
target_link_libraries(your_target PRIVATE adrius::adrius)
```

## Design

- **Zero-overhead backend abstraction** — algorithms are templates constrained by `adrius::Backend` (C++20 concept). No virtual dispatch; the compiler sees the full concrete type.
- **Two-stage API** — every algorithm exposes `preprocess_foo()` → `PreparedFoo` then `foo(PreparedFoo, params)`, plus a convenience overload that chains both. Lets you inspect, cache, or share the preprocessed lattice.
- **Named result structs** — no `std::tuple` with enum indices. Every return value has descriptive field names.
- **Value semantics** — results are returned by value (NRVO); `PreparedFoo` types are move-only. No owning pointers in the core API.
- **Lazy sequences** — `CFExpansionView` is a C++20 input range; composes with `std::views::take` etc.

See [`docs/design-decisions.md`](docs/design-decisions.md) for the full rationale behind every architectural choice.

## Documentation

| Document | What it covers |
|---|---|
| [`docs/getting-started.md`](docs/getting-started.md) | Installation, algorithm selection, params structs, backends, error handling |
| [`docs/lll.md`](docs/lll.md) | LLL lattice reduction — math, API, Lovász updates |
| [`docs/continued-fractions.md`](docs/continued-fractions.md) | CF expansion, convergents, best rational approximation |
| [`docs/simultaneous.md`](docs/simultaneous.md) | One-shot simultaneous Diophantine approximation |
| [`docs/illl-algorithm.md`](docs/illl-algorithm.md) | Iterated LLL (Bosma & Smeets 2010) — math, warm-start, guarantees |
| [`docs/boost-multiprecision-backend.md`](docs/boost-multiprecision-backend.md) | Arbitrary-precision backend — setup, precision options, performance |
| [`docs/design-decisions.md`](docs/design-decisions.md) | Full architectural rationale for every non-obvious choice |

## Dependencies

| Dependency | Version | How |
|---|---|---|
| C++20 compiler (MSVC 19.29+, GCC 11+, Clang 13+) | — | required |
| CMake | 3.21+ | required |
| Eigen3 | 3.4.0 | auto-fetched (or `ADRIUS_USE_SYSTEM_EIGEN=ON`) |
| GoogleTest | 1.14.0 | auto-fetched when `ADRIUS_BUILD_TESTS=ON` |

## License

MIT — see [LICENSE](LICENSE).  
Copyright (c) 2025 InGenifold Research LLC.
