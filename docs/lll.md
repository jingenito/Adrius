# LLL Lattice Basis Reduction

The Lenstra-Lenstra-Lovász (LLL) algorithm finds a "short, nearly orthogonal"
basis for an integer lattice.  It is the core primitive that ILLL and
`simultaneous_approx` build on.

## Mathematical background

Given an n×n real matrix B whose columns b₁, …, bₙ form a basis for a
lattice Λ = {∑ mᵢbᵢ : mᵢ ∈ ℤ}, LLL produces a new basis B' = B·U (where
U is a unimodular integer matrix) such that:

**Size-reduced**: for all i > j, |μᵢⱼ| ≤ η (default η = 0.51)

**Lovász condition**: for k = 2, …, n:
```
||b*ₖ||² ≥ (δ − μₖ,ₖ₋₁²) · ||b*ₖ₋₁||²
```

Here b*ₖ are the Gram-Schmidt orthogonalized vectors, μᵢⱼ = ⟨bᵢ, b*ⱼ⟩/||b*ⱼ||²
are the GSO coefficients, and δ ∈ (0.25, 1) is the Lovász parameter (default
0.75).

### What this guarantees

- The first basis vector b₁' satisfies ||b₁'|| ≤ 2^{(n-1)/4} · λ₁(Λ), where
  λ₁(Λ) is the length of the shortest non-zero lattice vector.
- Successive vectors satisfy ||b'ₖ|| ≤ 2^{(n-1)/2} · λₖ(Λ).
- The basis is output in O(n⁴ log B) arithmetic operations, where B bounds
  the input entries.

### Incremental GSO updates (Lovász)

Naively recomputing the full Gram-Schmidt decomposition after each column swap
costs O(n²) per swap and O(n⁵) total.  This implementation uses the Lovász
incremental update formulas (Cohen, §2.6.3) that update only the affected GSO
coefficients in O(n) per swap, giving the standard O(n⁴ log B) complexity.

---

## API

```cpp
#include <adrius/linalg/lll.hpp>

// Stage 1: validate the basis and compute the initial Gram-Schmidt decomposition.
// Throws DomainError if the basis is linearly dependent.
template <Backend B = DefaultBackend>
[[nodiscard]] PreparedLLL<B>
preprocess_lll(const matrix_of<B>& basis, LLLParams params = {});

// Stage 2: run LLL on a pre-computed state.
// Throws ConvergenceError if params.max_iter swaps are exceeded.
template <Backend B = DefaultBackend>
[[nodiscard]] LLLResult<B>
lll_reduce(PreparedLLL<B> prepared, LLLParams params = {});

// Convenience overload: chains both stages.
template <Backend B = DefaultBackend>
[[nodiscard]] LLLResult<B>
lll_reduce(const matrix_of<B>& basis, LLLParams params = {});
```

### Result types

```cpp
template <Backend B>
struct LLLResult {
    matrix_of<B>     reduced_basis;   // the LLL-reduced basis (columns)
    int_matrix_of<B> transform;       // unimodular U: reduced_basis = original · U
    std::size_t      swap_count;      // number of Lovász swaps (for profiling)
};

template <Backend B>
struct PreparedLLL {  // move-only
    matrix_of<B>     basis;
    int_matrix_of<B> transform;  // starts as identity
    GSOResult<B>     gso;
};
```

### Parameters

```cpp
struct LLLParams {
    double delta        = 0.75;   // Lovász δ ∈ (0.25, 1.0); closer to 1 → more reduced
    double eta          = 0.51;   // size-reduction threshold; 0.5 is the theoretical min
    std::size_t max_iter = 0;     // swap cap; 0 = unlimited
    GSOParams gso{};              // zero_threshold for linear-dependence detection
};
```

The classic LLL paper uses δ = 3/4 = 0.75, which is the default.  Values
closer to 1 yield shorter basis vectors at the cost of more swaps.

---

## Usage

### One-call reduction

```cpp
#include <adrius/linalg/lll.hpp>

adrius::EigenBackend::matrix_type basis(3, 3);
basis << 1, -1,  3,
         1,  0,  5,
         1,  2,  6;

auto result = adrius::lll_reduce<adrius::EigenBackend>(basis);

// result.reduced_basis — the LLL-reduced columns
// result.transform     — U such that reduced_basis = basis · U
// result.swap_count    — Lovász swaps performed
```

### Two-stage (inspect prepared state)

```cpp
auto prepared = adrius::preprocess_lll<adrius::EigenBackend>(basis);
// prepared.gso.norms_sq — initial squared GSO norms (diagnostic)

auto result = adrius::lll_reduce<adrius::EigenBackend>(std::move(prepared));
```

### Verifying the transform relationship

```cpp
// B' = B · U must hold exactly (up to floating-point).
double err = (basis * result.transform.cast<double>() - result.reduced_basis).norm();
// err should be < 1e-12 for well-conditioned inputs.
```

### Tuning the Lovász parameter

```cpp
// δ = 0.99: very high quality reduction, more swaps.
adrius::LLLParams params{ .delta = 0.99 };
auto result = adrius::lll_reduce<adrius::EigenBackend>(basis, params);
```

---

## Basis convention

Basis vectors are the **columns** of the matrix (Eigen column-major storage).
The unimodular transform acts on the **right**: `B' = B · U`.

This matches the Bosma-Smeets paper notation and LLL (1982).

---

## References

**Lenstra, A.K., Lenstra, H.W., Lovász, L.** (1982). *Factoring polynomials
with rational coefficients.* Mathematische Annalen, 261(4), 515–534.

**Cohen, H.** (1993). *A Course in Computational Algebraic Number Theory.*
Graduate Texts in Mathematics 138, Springer.  Algorithm 2.6.3 (size
reduction + Lovász swap), Algorithm 2.6.7 (incremental GSO update).

**Lovász, L.** (1986). *An Algorithmic Theory of Numbers, Graphs and
Convexity.* SIAM.  §1.3: incremental Gram-Schmidt updates after column swaps.
