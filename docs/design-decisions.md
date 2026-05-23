# Adrius — Design Decisions

This document records every non-obvious architectural choice made in this
library: what was decided, why, and what the alternative was. It exists so
that future contributors (and future Claude sessions) can make consistent
decisions in edge cases rather than guessing.

---

## 1. Backend abstraction: compile-time concepts, not runtime polymorphism

**Decision:** Backends are satisfied via C++20 `concept` constraints on a
struct. Algorithms are function templates parameterized on `Backend B`.

**Alternative rejected:** A pure-virtual `IBackend` interface with
`virtual` methods and `unique_ptr<IBackend>` ownership. This is the
"obvious OOP" answer and the wrong one here. Every algorithm call would
pay a virtual dispatch tax, and the compiler cannot inline or vectorize
across virtual boundaries — fatal for numerical code.

**Consequence:** Zero runtime overhead for the abstraction layer. The
compiler sees the full concrete type and can optimize freely. A future
runtime-swappable backend (e.g., for FFI/plugin use) would be a separate
`DynamicBackend` wrapper that type-erases internally; that is an explicit
future extension point, not the default path.

---

## 2. Result types: named-field structs, never tuples

**Decision:** Every algorithm returns a dedicated result struct with
descriptive field names. Example: `GSOResult<B>{ .Q = ..., .mu = ... }`
instead of `std::tuple<MatrixXd, MatrixXd>` with enum indices.

**Alternative rejected:** `std::tuple` + `enum class GSOType { Q, Mu }` for
index access. This pattern forces callers to remember positional semantics,
breaks autocomplete, and makes reading call sites nearly impossible without
the enum definition open.

**Convention:** All result structs carry `[[nodiscard]]` on the function
that produces them. All result structs have explicitly defaulted move
constructors (`= default`) so the compiler can NRVO or move freely.

---

## 3. Two-stage preprocessing API

**Decision:** Each algorithm exposes:
- `preprocess_foo(input, params)` → `PreparedFoo<B>` (validates, normalizes,
  constructs the lattice / GSO seed)
- `foo(PreparedFoo<B>, params)` → `FooResult<B>` (pure algorithm, no I/O)
- `foo(raw_input, params)` → `FooResult<B>` (convenience overload; chains
  the two stages internally)

**Why:** The old code fused input construction with algorithm execution in
functions like `SameDivisorFromRealVector()`. This made it impossible to:
- Test preprocessing in isolation from the algorithm.
- Inspect the constructed lattice before committing to a full run.
- Reuse a preprocessed state across multiple parameter sweeps.

**`PreparedXxx<B>` ownership model:** Move-only value type. Holds large
matrix data; copying it is expensive and almost never what the caller wants.
If a prepared state must be shared across multiple invocations on different
threads, wrap it in `std::shared_ptr<const PreparedFoo<B>>` at the call site.
The library does not impose shared ownership internally.

---

## 4. Memory model: value types and move semantics over smart pointers

**Decision:** Algorithm results and intermediate state are value types
returned by value. The compiler eliminates copies via NRVO. Smart pointers
are not used in the core API.

**Rationale:** For a header-only linear algebra library:
- Eigen matrices are already RAII types that manage their own heap storage.
- Returning `ILLLResult<B>` by value costs zero extra allocation when NRVO
  applies (single return path, no branches after construction).
- `unique_ptr<ILLLResult<B>>` adds an indirection and a heap allocation
  with no benefit.

**Where `std::shared_ptr<const T>` is appropriate:** When a caller computes
a `PreparedXxx<B>` once and wants to share it across multiple algorithm
invocations (e.g., a parameter sweep). This is a caller-side decision; the
library never forces it.

**Where `std::weak_ptr` is appropriate:** Paired with the above for cache
invalidation if a caching layer is ever added. Not in scope for v0.

**Non-owning views:** Function parameters that need a contiguous sequence
but do not take ownership use `std::span<const T>` (C++20). This accepts
`std::vector`, raw arrays, and `Eigen::Map` views without a copy.

---

## 5. Basis vector orientation: columns

**Decision:** Basis vectors are the **columns** of a `matrix_type`. The
`i`-th basis vector is `B::col(M, i)`.

**Why:** Eigen's default storage is column-major. The LLL and ILLL
literature (including Bosma & Smeets 2010) writes bases as column matrices.
Row-major orientation would require transposing every matrix operation,
fighting Eigen at every step.

**Consequence:** The GSO decomposition reads `M = Q * R` where columns of
`Q` are the orthogonalized vectors. The unimodular transform `U` acts on
the right: `M_new = M * U`.

---

## 6. Error handling: exceptions with a typed hierarchy

**Decision:** Algorithms signal failure by throwing from the `adrius::`
exception hierarchy (deriving `std::runtime_error`). No error-code returns,
no `std::optional` used to signal algorithmic failure.

**Alternative:** `tl::expected<T, E>` or a custom `Expected<T, E>`. This
would be the right call in C++23 (`std::expected` is standard), but adding a
third-party `expected` implementation as a dependency is not worth it for v0.
Revisit when the minimum standard is raised to C++23.

**Exception hierarchy (planned):**
```
adrius::Error              : std::runtime_error
  adrius::DomainError      — bad input (e.g., linearly dependent basis)
  adrius::ConvergenceError — algorithm failed to converge within bounds
  adrius::PrecisionError   — floating-point precision insufficient for input
```

**Note:** All `noexcept` annotations must be accurate. Internal helpers that
cannot throw are marked `noexcept`. Algorithm entry points are not `noexcept`.

---

## 7. Backend type system: single Backend carries both scalar and integer types

**Decision:** The `Backend` concept requires both `scalar_type` (floating-
point, e.g. `double`) and `integer_type` (signed integer, e.g. `int64_t`).
A single template parameter `B` controls both.

**Alternative rejected:** Two separate Backend template parameters
`<Backend Float, Backend Int>`. This is overcomplicated for the actual use
case: in ILLL, the floating-point GSO coefficients and the integer unimodular
transform are always paired. Splitting them into independent type parameters
would make every algorithm signature noisier with no real benefit.

**Integer precision policy (v0):** `integer_type` defaults to `int64_t`.
Overflow is theoretically possible on inputs with very large partial quotients
or many ILLL iterations. A `MpBackend` using `__int128` or Boost.Multiprecision
as `integer_type` is a natural future extension; the concept is already
designed to accommodate it.

---

## 8. Public namespace: flat `adrius::`

**Decision:** All public entry points live in the `adrius::` namespace. No
`adrius::approx::`, `adrius::linalg::`, etc. in the public API.

**Alternative rejected:** Nested namespaces like
`adrius::approx::simultaneous::illl()`. The old code had
`SimultaneousDiophantine::IteratedLLL_Dyadic()` — verbose, hard to type,
confusing in autocomplete.

**Internal namespaces:** `adrius::detail::` for implementation helpers not
intended for public use. These are not subject to API stability guarantees.

---

## 9. Parameter structs for all tuneable algorithms

**Decision:** Any algorithm with more than one tuning parameter takes a
`FooParams` struct with named fields and default values. No long positional
parameter lists.

**Example:**
```cpp
struct LLLParams {
    double delta        = 0.75;  // Lovász condition threshold (0.25, 1.0)
    double eta          = 0.51;  // size-reduction threshold
    std::size_t max_iter = 0;    // 0 = unlimited
};
```

Callers can use designated initializers: `LLLParams{ .delta = 0.99 }`.

---

## 10. Incremental GSO updates in LLL (Lovász update formulas)

**Decision:** The LLL implementation uses incremental GSO updates after each
column swap, not a full recomputation.

**Why:** A naive LLL recomputes the full Gram-Schmidt orthogonalization after
every swap — O(n²) per swap, O(n³) total swaps in the worst case, giving
O(n⁵) overall. The Lovász update formulas reduce each swap's GSO update to
O(n), giving the standard O(n⁴ log B) complexity.

**Reference:** Lovász (1986), "An Algorithmic Theory of Numbers, Graphs and
Convexity", Section 1.3. See also Cohen, "A Course in Computational Algebraic
Number Theory", Algorithm 2.6.3.

---

## 11. Lazy evaluation for continued fraction sequences

**Decision:** `cf_expansion()` returns a lazy range view
(`CFExpansionView`) rather than a pre-allocated `std::vector<int64_t>`.

**Why:** The caller decides how many partial quotients they need. A lazy view:
- Computes quotients on demand; no allocation until materialized.
- Composes with `std::ranges` algorithms (`| std::views::take(N)`).
- Allows early termination (e.g., stop when the convergent is within bound)
  without wasted computation.

**Materialization:** `adrius::collect(view)` or `std::ranges::to<std::vector>`
(C++23) converts to a concrete container when needed.

---

## 12. `Eigen::Map<>` for zero-copy user data ingestion

**Decision:** The EigenBackend exposes `map_vector()` and `map_matrix()`
factory functions that wrap user-provided `double[]` data using
`Eigen::Map<>` without copying.

**Why:** A researcher who already has data in a C array or a contiguous
buffer should not pay for an allocation and copy just to call an algorithm.
`Eigen::Map` is a zero-overhead wrapper that presents raw memory as an
Eigen type, satisfying the backend's `vector_type` requirements.

**Lifetime contract:** The `Map` is valid only as long as the underlying
buffer lives. This is documented on the factory function and is the caller's
responsibility.

---

## 13. `[[nodiscard]]` policy

**Decision:** Applied to:
- All functions returning result structs (`GSOResult`, `LLLResult`, etc.)
- All `preprocess_*` functions
- `Rational<Int>` arithmetic operators

**Not applied to:** void functions, functions called for side effects.

---

## 14. Install / `find_package` support from day one

**Decision:** The CMake build exports proper config files so
`find_package(adrius 0.1 REQUIRED)` works for system-installed consumers.
FetchContent is also supported.

**Why:** Retrofitting install support onto a CMake project that wasn't
designed for it is painful. The `$<BUILD_INTERFACE:...>/$<INSTALL_INTERFACE:...>`
generator expressions in `target_include_directories` are already in place;
adding `AdriusConfig.cmake.in` costs ~30 lines now vs a painful refactor later.

**`ADRIUS_ENABLE_INSTALL` guard:** Defaults to `ON` only when
`PROJECT_IS_TOP_LEVEL` is true. Prevents polluting the install prefix of
parent projects that consume Adrius via FetchContent.

---

## 15. Copyright header on every file

Every source file begins with:
```cpp
// Copyright (c) 2025 InGenifold Research LLC. MIT License.
```

No exceptions. This is a CI-enforceable policy (grep-based check in the
future).
