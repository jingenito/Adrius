# Simultaneous Diophantine Approximation

`simultaneous_approx` finds a single common denominator q and integer
numerators p₁, …, pₙ such that |q·αᵢ − pᵢ| is small for all i
simultaneously.  This is the one-shot variant; for an iterated sequence with
improving quality use `illl` (see `docs/illl-algorithm.md`).

## Mathematical background

**Dirichlet's simultaneous approximation theorem**: For any α = (α₁, …, αₙ)
and any Q ≥ 1, there exists an integer q with 1 ≤ q ≤ Qⁿ and integers
p₁, …, pₙ such that
```
|q·αᵢ − pᵢ| ≤ 1/Q   for all i = 1, …, n.
```

The LLL-based approach finds such a (q, p) efficiently by embedding the
problem in a lattice and reducing it.

### Lattice construction

Given α = (α₁, …, αₙ) and scale N, build the (n+1)×(n+1) matrix:
```
L = [ N    0   0  … 0 ]
    [ Nα₁  1   0  … 0 ]
    [ Nα₂  0   1  … 0 ]
    [  ⋮         ⋱   ⋮ ]
    [ Nαₙ  0   0  … 1 ]
```

The first column of the LLL-reduced basis of L has the form
(q·N, p₁, …, pₙ) for integers q, pᵢ.  Dividing by N extracts q; the
remaining entries are the numerators pᵢ.

By Dirichlet, this produces approximation quality ≈ 1/N^{1/n}.  Larger N
gives smaller errors but may require more LLL work.

### Relation to ILLL

`simultaneous_approx` performs **one LLL call** at a fixed scale N and
returns a single (q, p).  `illl` iterates: it repeatedly applies LLL with a
geometrically shrinking scale, producing a sequence of (q, p) pairs with
improving quality and growing denominators, with the Bosma-Smeets bound on the
Dirichlet coefficient q^{1/n}·max|q·αᵢ−pᵢ|.

Use `simultaneous_approx` when you want one approximation at a specific scale.
Use `illl` when you want a sequence of improving approximations.

---

## API

```cpp
#include <adrius/approx/simultaneous.hpp>

// Stage 1: build the lattice for α at scale N.
// Throws DomainError if alpha is empty or N ≤ 0.
template <Backend B = DefaultBackend>
[[nodiscard]] PreparedSimultaneous<B>
preprocess_simultaneous(std::span<const scalar_of<B>> alpha,
                        scalar_of<B> N);

// Stage 2: LLL-reduce and extract (q, p₁, …, pₙ).
template <Backend B = DefaultBackend>
[[nodiscard]] SimultApproxResult<B>
simultaneous_approx(PreparedSimultaneous<B> prepared, LLLParams params = {});

// Convenience overload: chains both stages.
template <Backend B = DefaultBackend>
[[nodiscard]] SimultApproxResult<B>
simultaneous_approx(std::span<const scalar_of<B>> alpha,
                    scalar_of<B> N,
                    LLLParams params = {});
```

### Result types

```cpp
template <Backend B>
struct SimultApproxResult {
    integer_of<B>              denominator;  // the common denominator q
    std::vector<integer_of<B>> numerators;   // [p₁, …, pₙ]
    scalar_of<B>               quality;      // max_i |q·αᵢ − pᵢ|
};

template <Backend B>
struct PreparedSimultaneous {  // move-only
    matrix_of<B>              lattice;
    std::vector<scalar_of<B>> alpha;   // original targets (for quality computation)
    scalar_of<B>              scale;   // N
};
```

---

## Usage

### One-call approximation

```cpp
#include <adrius/approx/simultaneous.hpp>
#include <numbers>
#include <span>

const std::vector<double> alpha = {std::numbers::sqrt2, std::numbers::pi};

auto result = adrius::simultaneous_approx<adrius::EigenBackend>(
    std::span<const double>{alpha},
    /*scale N=*/ 1e6);

// result.denominator — q
// result.numerators  — {p₁, p₂}
// result.quality     — max|q·αᵢ − pᵢ|
std::cout << "q = "    << result.denominator    << '\n'
          << "p₁ = "   << result.numerators[0]  << '\n'
          << "p₂ = "   << result.numerators[1]  << '\n'
          << "err = "  << result.quality         << '\n';
```

### Choosing the scale N

The expected quality is ≈ 1/N^{1/n}.  For n=2 and quality target 10⁻⁴,
use N ≈ (10⁴)² = 10⁸.

```cpp
// For n = 2 targets, quality ≈ 1/N^{0.5}:
// N = 1e4 → quality ≈ 0.01
// N = 1e6 → quality ≈ 0.001
// N = 1e8 → quality ≈ 0.0001
```

### Two-stage (inspect the lattice before reduction)

```cpp
auto prepared = adrius::preprocess_simultaneous<adrius::EigenBackend>(
    std::span<const double>{alpha}, /*N=*/ 1e6);

// Inspect constructed lattice:
std::cout << "Lattice:\n";
for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j)
        std::cout << adrius::EigenBackend::get(prepared.lattice, i, j) << ' ';
    std::cout << '\n';
}

// Run LLL with custom params:
adrius::LLLParams lll_params{ .delta = 0.99 };
auto result = adrius::simultaneous_approx<adrius::EigenBackend>(
    std::move(prepared), lll_params);
```

### Parameter sweep over different scales

```cpp
for (double log_N : {4.0, 5.0, 6.0, 7.0}) {
    const double N = std::pow(10.0, log_N);
    auto result = adrius::simultaneous_approx<adrius::EigenBackend>(
        std::span<const double>{alpha}, N);
    std::cout << "N=1e" << log_N
              << "  q=" << result.denominator
              << "  err=" << result.quality << '\n';
}
```

---

## Comparison with ILLL

| Property | `simultaneous_approx` | `illl` |
|---|---|---|
| Number of LLL calls | 1 | many (iterates) |
| Quality control | fixed by N | improves per iteration |
| Output | one (q, p) | sequence of (q, p) with growing q |
| Quality bound | ~1/N^{1/n} (Dirichlet) | q^{1/n}·quality ≤ C (Bosma-Smeets) |
| Parameters | scale N | epsilon, max_iterations, max_denominator |

---

## References

**Dirichlet, P.G.L.** (1842). Verallgemeinerung eines Satzes aus der Lehre
von den Kettenbrüchen nebst einigen Anwendungen auf die Theorie der Zahlen.
*Bericht Königl. Preuss. Akad. Wiss.*, 93–95.

**Lagarias, J.C.** (1985). *The computational complexity of simultaneous
Diophantine approximation problems.* SIAM Journal on Computing, 14(1), 196–209.
LLL-based simultaneous approximation algorithm and complexity bounds.

**Bosma, W. & Smeets, I.** (2010). *Finding simultaneous Diophantine
approximations with prescribed quality.* arXiv:1001.4455.  The ILLL algorithm
that extends this one-shot approach to an iterative sequence.
