# Continued Fractions and Best Rational Approximation

This document covers `cf_view`, `cf_expansion`, `cf_convergents`, and
`best_rational` — the four functions in `adrius/approx/continued_fraction.hpp`
and `adrius/approx/rational.hpp`.

## Mathematical background

Every real number x has a **continued fraction expansion**:
```
x = a₀ + 1/(a₁ + 1/(a₂ + 1/(a₃ + ⋯)))  =:  [a₀; a₁, a₂, a₃, …]
```

The integers aₖ ≥ 1 (a₀ may be any integer) are the **partial quotients**.
For irrational x the sequence is infinite; for rational x it terminates.

### Convergents

Truncating the CF at depth k gives the **k-th convergent** pₖ/qₖ, computed
by the three-term recurrences:
```
p₋₁ = 1,   p₀ = a₀,   pₖ = aₖ·pₖ₋₁ + pₖ₋₂
q₋₁ = 0,   q₀ = 1,    qₖ = aₖ·qₖ₋₁ + qₖ₋₂
```

Convergents are the **best rational approximations** to x: no rational p/q
with q ≤ qₖ approximates x better than pₖ/qₖ.

### Dirichlet's theorem

For any irrational x and any Q ≥ 1, there exists p/q with 1 ≤ q ≤ Q such
that |x − p/q| < 1/qQ.  Convergents achieve this bound; `best_rational`
finds exactly such a fraction.

---

## API

```cpp
#include <adrius/approx/continued_fraction.hpp>
#include <adrius/approx/rational.hpp>
```

### Lazy range: `cf_view`

```cpp
CFExpansionView cf_view(double x, CFParams params = {});
```

Returns a C++20 input range that yields partial quotients on demand.
Composes with `std::views::` adaptors.

```cpp
// Print the first 8 partial quotients of π = [3; 7, 15, 1, 292, 1, 1, 1, …]
for (std::int64_t a : adrius::cf_view(std::numbers::pi) | std::views::take(8))
    std::cout << a << ' ';   // 3 7 15 1 292 1 1 1
```

### Eager vector: `cf_expansion`

```cpp
std::vector<std::int64_t> cf_expansion(double x, CFParams params = {});
```

Materializes the entire expansion up to `params.max_depth` quotients.

```cpp
auto qs = adrius::cf_expansion(std::numbers::phi, {.max_depth = 12});
// qs = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}  (φ = [1; 1, 1, 1, …])
```

### Convergents: `cf_convergents`

```cpp
std::vector<Convergent> cf_convergents(double x, CFParams params = {});

struct Convergent {
    std::int64_t p;  // numerator
    std::int64_t q;  // denominator (always > 0)
};
```

Returns all convergents up to `params.max_depth`.

```cpp
auto convs = adrius::cf_convergents(std::numbers::pi, {.max_depth = 5});
// convs = [{3,1}, {22,7}, {333,106}, {355,113}, {103993,33102}]
```

### Best rational: `best_rational`

```cpp
Rational<std::int64_t> best_rational(double x, RationalApproxParams params = {});

struct RationalApproxParams {
    std::int64_t max_denominator = 1'000'000LL;
};
```

Returns the rational p/q with q ≤ `max_denominator` that minimizes |x − p/q|.
Uses convergents and Farey semi-convergents (Hardy & Wright, Theorem 171).

```cpp
auto r = adrius::best_rational(std::numbers::pi, {.max_denominator = 1000});
// r = 355/113 — within 2.7e-7 of π, the best possible with q ≤ 1000
std::cout << r;                 // 355/113
std::cout << r.numerator();     // 355
std::cout << r.denominator();   // 113
std::cout << r.to_double();     // 3.14159292...
```

---

## The `Rational<Int>` type

`best_rational` returns `Rational<std::int64_t>`, defined in
`adrius/util/rational_type.hpp`.  It is always in reduced form (gcd = 1,
denominator > 0).

```cpp
adrius::Rational<std::int64_t> r{22, 7};

r.numerator();    // 22
r.denominator();  // 7
r.to_double();    // 3.142857...
r.is_integer();   // false

// Arithmetic (returns reduced form):
auto s = r + adrius::Rational<std::int64_t>{1, 7};  // 23/7
auto t = r * r;                                      // 484/49
```

---

## Parameters

```cpp
struct CFParams {
    std::size_t max_depth       = 64;    // maximum number of partial quotients
    double      zero_threshold  = 1e-14; // fractional parts below this → exact rational
};
```

Decrease `max_depth` to limit computation; `max_depth = 90` is sufficient to
reach any `int64_t` denominator (convergent denominators grow at least as fast
as Fibonacci numbers).

---

## Common patterns

### Early termination with `std::views::take_while`

```cpp
// Stop when the partial quotient exceeds 100 (large quotients signal near-rationals).
for (auto a : adrius::cf_view(x) | std::views::take_while([](auto q){ return q <= 100; }))
    process(a);
```

### Finding the convergent index where q first exceeds a bound

```cpp
auto convs = adrius::cf_convergents(x);
auto it = std::ranges::find_if(convs, [](const auto& c){ return c.q > 1000; });
auto best_before = *std::prev(it);  // last convergent with q ≤ 1000
```

### Checking approximation quality

```cpp
for (const auto& [p, q] : adrius::cf_convergents(std::numbers::e, {.max_depth = 8})) {
    double err = std::abs(std::numbers::e - static_cast<double>(p)/q);
    // err < 1/(q · q_next) for every convergent
}
```

---

## References

**Hardy, G.H. & Wright, E.M.** (1979). *An Introduction to the Theory of
Numbers*, 5th ed. Oxford University Press.  Chapter X (continued fractions),
Theorem 171 (best approximation property of convergents), Theorem 181 (Hurwitz
bound q·|x − p/q| < 1/√5 infinitely often).

**Khinchin, A.Ya.** (1964). *Continued Fractions.* University of Chicago
Press.  §1–2: three-term recurrences and convergent error bounds.
