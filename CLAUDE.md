# Adrius — Diophantine Approximation Library
**InGenifold Research LLC | MIT License**

## Project Overview
C++20 header-only library for Diophantine approximation (ILLL, LLL, continued
fractions, simultaneous approximation). Pluggable linear algebra backends;
Eigen3 is the default.

## Directory Layout
- `include/adrius/` — all public headers (header-only)
  - `core/`    — concepts, traits, params structs (zero external dependencies)
  - `backend/` — backend implementations; only file allowed to include Eigen
  - `linalg/`  — Gram-Schmidt, LLL (depend only on `core/`)
  - `approx/`  — continued fractions, rational, simultaneous, ILLL
  - `util/`    — Rational<Int> and other standalone utilities
- `tests/`   — GoogleTest suite; one file per algorithm
- `examples/` — standalone consumer programs
- `cmake/`   — `AdriusConfig.cmake.in`, version file
- `docs/`    — `design-decisions.md` (authoritative record of architectural choices)

## Code Style Rules
- C++20 standard; no C++23 features until the minimum is explicitly raised
- NEVER use `using namespace` at file scope; always qualify or use block-scope `using`
- No raw owning pointers. Prefer value types and move semantics over smart pointers;
  see "Memory model" below for when smart pointers are appropriate
- Template code lives in headers; no explicit instantiation unless justified
- `[[nodiscard]]` on every function that returns a result struct or prepared state
- Named-field structs for all multi-value returns — never `std::tuple` with enum indices
- Every algorithm with more than one tuning parameter takes a `FooParams` struct with
  named fields and sensible defaults; callers use designated initializers
- `adrius::detail::` namespace for internal helpers; never expose in public headers
- Copyright header on every file (first line, before any includes):
  // Copyright (c) 2025 InGenifold Research LLC. MIT License.

## Backend Abstraction Rules
- Algorithm templates must only `#include` headers from `core/`. They must never
  directly include Eigen or any other backend header.
- `backend/eigen.hpp` is the ONLY file permitted to `#include <Eigen/...>`.
- Every backend implementation must satisfy `static_assert(adrius::Backend<MyBackend>)`
  to catch regressions at compile time.
- Basis vectors are the **columns** of `matrix_type` (matches Eigen column-major storage
  and Bosma-Smeets paper notation). The unimodular transform acts on the right: `M' = M·U`.

## Memory Model
Prefer value types returned by value (NRVO eliminates copies). Smart pointer guidance:
- `unique_ptr`: use only for owning polymorphic types or Pimpl — not applicable in the
  current header-only, concept-based design. Reserve for a future runtime-swappable backend.
- `shared_ptr<const PreparedFoo<B>>`: appropriate when a caller shares a preprocessed state
  across multiple algorithm invocations (e.g., a parameter sweep). The library never
  imposes shared ownership internally.
- `weak_ptr`: paired with the above for cache invalidation if a caching layer is added.
- `std::span<const T>`: use for non-owning sequence parameters; accepts vectors, arrays,
  and Eigen maps without a copy.
- `Eigen::Map<>`: zero-copy wrapping of user-provided raw buffers (caller owns lifetime).

## API Design Rules
- Public entry points live in flat `adrius::` namespace only. No `adrius::approx::illl()`.
- Every algorithm exposes a two-stage API:
    1. `preprocess_foo(input, params)` → `PreparedFoo<B>` (move-only value type)
    2. `foo(PreparedFoo<B>, params)` → `FooResult<B>`
  plus a convenience overload `foo(raw_input, params)` that chains both stages.
- `PreparedFoo<B>` types are move-only. If sharing across invocations, the caller wraps
  in `shared_ptr<const PreparedFoo<B>>`.
- Algorithms signal failure by throwing from the `adrius::` exception hierarchy:
    adrius::Error : std::runtime_error
      adrius::DomainError      — bad input (linearly dependent basis, out-of-range param)
      adrius::ConvergenceError — algorithm did not converge within specified bounds
      adrius::PrecisionError   — floating-point precision insufficient for the input
- `integer_type` defaults to `int64_t`. No arbitrary-precision integers in v0.

## Algorithm Implementation Rules
- GSO updates inside LLL must use Lovász incremental update formulas — never a full
  recompute after a column swap. (See `docs/design-decisions.md` §10 for reference.)
- Continued fraction sequences are lazy range views (`CFExpansionView`), not
  pre-allocated `std::vector`. Materialize only when the caller requests it.
- Every algorithm must have at minimum one test verifying a hand-computed small example
  before any randomized/property-based tests are added.

## Build System
CMake 3.21+ with FetchContent for all dependencies.
- `cmake -B build -DCMAKE_BUILD_TYPE=Release`
- `cmake --build build`
- Run tests: `ctest --test-dir build --output-on-failure`
- Dependencies auto-fetched: Eigen3, GoogleTest
- `ADRIUS_USE_SYSTEM_EIGEN=ON` to use an installed Eigen instead of fetching
- `ADRIUS_ENABLE_INSTALL` guards install targets; defaults ON only when top-level project
- `find_package(adrius)` is supported via `cmake/AdriusConfig.cmake.in`

## Key Algorithms
Primary target: ILLL from Bosma & Smeets (2010).
Previous scratch implementation had bugs in the bounds at Lemma 25 (eq. 28)
and Theorem 23 (eq. 24) — do not port any code from the old MUSE repository.
The old repo (github.com/jingenito/MUSE) is reference only for consumer-facing
API patterns (CPP, CPPTestBench, ProjectBosmaSmeets), not for implementation.

## Platform Targets
Linux (GCC/Clang), Windows (MSVC), macOS (AppleClang).
Use CMake abstractions; no platform-specific `#ifdef`s unless truly unavoidable.
MinGW is a secondary target; do not break it but do not optimize for it.

## Git Workflow & Branching

**Always create a new branch off the latest `main` when starting work:**

```bash
git checkout main
git pull origin main
git checkout -b <branch-name>
```

### Branch Naming
- `feature/description` — new feature (e.g., `feature/simultaneous-approx`)
- `fix/description` — bug fix (e.g., `fix/lll-convergence-edge-case`)
- `docs/description` — documentation/guides (e.g., `docs/instruction-files`)
- `refactor/description` — code cleanup, no behavior change (e.g., `refactor/backend-abstraction`)
- `test/description` — tests and test infrastructure (e.g., `test/property-based-lll`)

### Commit Messages
- **First line**: Imperative mood, < 70 chars, starts with verb
  - ✅ "Implement ILLL algorithm from Bosma & Smeets (2010)"
  - ✅ "Add comprehensive ILLL algorithm documentation"
  - ❌ "ILLL implementation" (not imperative)
  - ❌ "This commit implements the ILLL algorithm..." (too wordy)
- **Body** (if needed): Detailed explanation of *why*, not *what*
  - Explain trade-offs, constraints, design decisions
  - Reference issue numbers, papers, or related commits
- **Sign-off**: Always include `Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>`

Example:
```
Implement ILLL algorithm from Bosma & Smeets (2010) with correct scale update

- Add ILLLResult<B> struct for storing approximation relations
- Implement preprocess_illl<B>() building Bosma-Smeets lattice with c₀=ε^{n+1}
- Implement illl<B>() iteration with scale_ratio = 2^{-n(n+1)/4}

Fixes issue with denominator extraction by using balanced lattice construction.

Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>
```

### Pull Request Workflow
1. **Before opening PR**: Run full test suite locally, verify all tests pass
2. **PR description**: Use `.github/pull_request_template.md` format
   - Summary of what changed and why
   - Type of change (Feature/Fix/Docs/Refactor)
   - Test results proof
   - Research papers referenced (if applicable)
3. **After merge**: Branch is deleted (locally and remotely)

## Documentation Requirements

### Algorithm or API Changes
- **Must include** `docs/algorithm-name.md` with:
  - Overview and purpose
  - Mathematical foundation (with paper references)
  - Implementation details and key decisions
  - Usage examples
  - References section citing all papers used
- **Must reference papers** in PR description with specific sections/equations used
- **Must update** README.md if adding new public API or changing behavior

### Code Comments
- No comments explaining *what* the code does (variable names and function names do that)
- Only comments explaining *why*:
  - Non-obvious design decisions
  - Subtle invariants or constraints
  - Workarounds for specific bugs
  - References to papers/algorithms (e.g., "Lovász incremental update from LLL paper §3")
- Example: ✅ "Rescale only column 0; keeps basis nearly reduced for amortized cost"
- Example: ❌ "Rescale column 0 by the scale ratio"

### Test Requirements
- Every algorithm must have at minimum one hand-computed example test
- Tests should verify both correctness and bounds
- Test names describe what is being verified, not just "Test1"
- Example: ✅ `TEST(ILLL, OneDimensionalBoundedDirichletCoefficient)`
- Example: ❌ `TEST(ILLL, Test1)`

## Agile Lifecycle Rules
1. **In Dev**: Always create a feature branch off the latest `main`. Implement changes and author unit/integration tests. See "Git Workflow & Branching" above.
2. **Self-Testing & QA Gating**: Run the full local test suite (`ctest --test-dir build --output-on-failure`). Resolve all compilation errors, linter issues, and test failures before submitting for review.
3. **Ready for Review**: Push branch and create PR. Include:
   - Summary of changes and why they were made
   - Type of change (Feature/Fix/Docs/Refactor)
   - Proof of passing test suite (ctest output or screenshot)
   - Research papers referenced (with specific sections if applicable)
   - For algorithm changes: link to `docs/algorithm-name.md`