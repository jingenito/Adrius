# Getting Started with Adrius

Adrius is a C++20 header-only library for Diophantine approximation.  This
guide walks through adding it to a project, choosing the right algorithm, and
understanding the patterns common to the entire API.

## Adding Adrius to a project

### Via FetchContent (recommended)

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

No other steps required.  Eigen3 is fetched automatically.

### After system install

```cmake
find_package(adrius 0.1 REQUIRED)
target_link_libraries(your_target PRIVATE adrius::adrius)
```

Install with `cmake --install build` from your build directory.

### Including headers

```cpp
// Everything with the default EigenBackend:
#include <adrius/adrius.hpp>

// Or include only what you need:
#include <adrius/approx/illl.hpp>
#include <adrius/linalg/lll.hpp>
```

---

## Choosing an algorithm

| Goal | Algorithm | Header |
|---|---|---|
| Approximate one real number p/q with q ≤ Q | `best_rational` | `approx/rational.hpp` |
| Continued-fraction expansion of one real | `cf_view` / `cf_expansion` | `approx/continued_fraction.hpp` |
| One-shot simultaneous approximation at fixed scale | `simultaneous_approx` | `approx/simultaneous.hpp` |
| Iterated simultaneous approximation with improving quality | `illl` | `approx/illl.hpp` |
| Find a short basis for a lattice | `lll_reduce` | `linalg/lll.hpp` |

As a rule of thumb:
- For **one** target, use `best_rational` or `cf_view`.
- For **multiple** targets with a fixed scale, use `simultaneous_approx`.
- For **multiple** targets where you want increasingly good approximations over
  many iterations (and the Bosma-Smeets quality bound), use `illl`.

---

## The two-stage API

Every non-trivial algorithm in Adrius exposes three entry points:

```cpp
// Stage 1: validate input and construct the lattice / initial state.
auto prepared = adrius::preprocess_foo<B>(input, params);

// Stage 2: run the algorithm on the prepared state.
auto result = adrius::foo<B>(std::move(prepared), params);

// Convenience: chain both stages in one call.
auto result = adrius::foo<B>(input, params);
```

`PreparedFoo<B>` is a **move-only** value type.  Copy it into a
`std::shared_ptr<const PreparedFoo<B>>` if you need to share it across
multiple calls (e.g., a parameter sweep over `params`).

---

## Params structs

All tunable parameters are passed via named-field structs with sensible
defaults.  Use C++20 designated initializers to override only what you need:

```cpp
adrius::ILLLParams params{
    .max_iterations  = 50,
    .max_denominator = 1'000'000LL,
    .quality_tol     = 1e-10,
};
```

The full parameter structs are defined in `adrius/core/params.hpp`:

```cpp
struct GSOParams   { double zero_threshold = 1e-14; };

struct LLLParams   { double delta = 0.75;
                     double eta   = 0.51;
                     std::size_t max_iter = 0;
                     GSOParams gso{}; };

struct ILLLParams  { LLLParams lll{};
                     std::size_t max_iterations  = 30;
                     std::int64_t max_denominator = 1'000'000LL;
                     double epsilon     = 0.5;
                     double quality_tol = 1e-12; };

struct CFParams    { std::size_t max_depth       = 64;
                     double zero_threshold = 1e-14; };

struct RationalApproxParams { std::int64_t max_denominator = 1'000'000LL; };
```

---

## Backends

The `B` template parameter selects the linear algebra backend.

### EigenBackend (default)

Uses `double` precision and Eigen3 matrices.  No extra dependencies.

```cpp
#include <adrius/adrius.hpp>

auto result = adrius::illl<adrius::EigenBackend>(alpha, params);
// Or equivalently, EigenBackend is the DefaultBackend:
auto result = adrius::illl(alpha, params);
```

`EigenBackend` type aliases:

| Alias | Concrete type |
|---|---|
| `EigenBackend::scalar_type` | `double` |
| `EigenBackend::integer_type` | `std::int64_t` |
| `EigenBackend::matrix_type` | `Eigen::MatrixXd` |
| `EigenBackend::vector_type` | `Eigen::VectorXd` |
| `EigenBackend::int_matrix_type` | `Eigen::MatrixXi` |

### BoostBackend (arbitrary precision)

Uses Boost.Multiprecision for ill-conditioned inputs or very large denominators.
Requires `-DADRIUS_ENABLE_BOOST_MULTIPRECISION=ON` at CMake configure time.

```cpp
#include <adrius/adrius.hpp>

using B = adrius::BoostBackend<50>;  // 50 decimal digits
using Scalar = B::scalar_type;       // boost::multiprecision::number<...>

std::vector<Scalar> alpha = {
    Scalar("1.41421356237309504880168872420969807856967187537694"),
    Scalar("3.14159265358979323846264338327950288419716939937510"),
};

auto result = adrius::illl<B>(std::span<const Scalar>{alpha},
                              adrius::ILLLParams{.max_denominator = 1'000'000'000LL});
```

See `docs/boost-multiprecision-backend.md` for precision options and
performance guidance.

### Custom backend

Any type satisfying `adrius::Backend` works.  See `docs/design-decisions.md`
§1 and the concept definition in `include/adrius/core/concepts.hpp`.

```cpp
static_assert(adrius::Backend<MyBackend>);  // fails loudly if incomplete
```

---

## Error handling

All errors throw from the `adrius::Error` hierarchy (derives
`std::runtime_error`):

```cpp
try {
    auto result = adrius::illl(alpha, params);
} catch (const adrius::DomainError& e) {
    // Bad input: empty alpha, epsilon out of (0,1), etc.
} catch (const adrius::ConvergenceError& e) {
    // LLL swap limit exceeded (params.lll.max_iter > 0 and was reached).
} catch (const adrius::Error& e) {
    // Any other Adrius error.
}
```

Normal algorithmic non-convergence (denominator cap reached, quality_tol met,
iteration limit) is **not** an exception — the algorithm returns whatever
relations it found and you inspect `result.iterations`.

---

## First example: best rational approximation

```cpp
#include <adrius/adrius.hpp>
#include <iostream>
#include <numbers>

int main() {
    // 355/113 — best approximation to π with denominator ≤ 1000.
    auto r = adrius::best_rational(std::numbers::pi, {.max_denominator = 1000});
    std::cout << r << '\n';          // 355/113
    std::cout << r.to_double() << '\n';
    std::cout << r.numerator() << " / " << r.denominator() << '\n';
    return 0;
}
```

## Second example: iterated simultaneous approximation

```cpp
#include <adrius/adrius.hpp>
#include <iostream>
#include <numbers>
#include <span>

int main() {
    const std::vector<double> alpha = {std::numbers::sqrt2, std::numbers::pi};

    adrius::ILLLParams params{
        .max_iterations  = 20,
        .max_denominator = 100'000,
        .quality_tol     = 1e-8,
    };

    auto result = adrius::illl<adrius::EigenBackend>(
        std::span<const double>{alpha}, params);

    for (std::size_t k = 0; k < result.iterations; ++k) {
        const auto& rel = result.relations[k];   // [q, p1, p2]
        std::cout << "q=" << rel[0]
                  << " p1=" << rel[1]
                  << " p2=" << rel[2]
                  << " quality=" << result.quality[k] << '\n';
    }
    return 0;
}
```

---

## Python bindings

All public algorithms are also available from Python via optional pybind11
bindings.  Build with `-DADRIUS_BUILD_PYTHON=ON` and see
`docs/python-bindings.md` for the full API reference, build instructions, and
test-running guide.

---

## Further reading

- `docs/lll.md` — LLL lattice basis reduction
- `docs/continued-fractions.md` — CF expansion and best rational approximation
- `docs/simultaneous.md` — One-shot simultaneous Diophantine approximation
- `docs/illl-algorithm.md` — Iterated LLL (Bosma & Smeets 2010)
- `docs/boost-multiprecision-backend.md` — Arbitrary-precision backend
- `docs/design-decisions.md` — Full architectural rationale
