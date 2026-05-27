# Python Bindings

Adrius ships optional Python bindings built with pybind11.  The bindings
expose all public algorithms via the `adrius` Python package and assume the
`EigenBackend` (double precision).  NumPy arrays are accepted wherever the C++
API takes a span of doubles.

---

## Build

The bindings are off by default and do not affect C++-only users.

```bash
cmake -B build -DADRIUS_BUILD_PYTHON=ON \
               -DADRIUS_BUILD_TESTS=OFF \
               -DADRIUS_BUILD_EXAMPLES=OFF
cmake --build build --target _adrius
```

The compiled extension (`_adrius.pyd` on Windows, `_adrius.so` on Linux/macOS)
is placed inside `build/python/adrius/` alongside the pure-Python `__init__.py`.

### Using without installing

Point `PYTHONPATH` at the build output:

```bash
export PYTHONPATH=build/python   # Linux / macOS
$env:PYTHONPATH = "build\python" # PowerShell
```

### pybind11 dependency

The build prefers a system or pip-installed pybind11.  If none is found it
fetches pybind11 v2.13.6 via FetchContent automatically.

```bash
pip install pybind11   # optional: speeds up configure step
```

---

## Quick start

```python
import adrius
import numpy as np

# Best rational approximation
r = adrius.best_rational(3.14159265, max_denominator=1000)
print(r)          # 355/113
print(float(r))   # 3.1415929203539825

# Continued-fraction expansion
print(adrius.cf_expansion(3.14159265358979, max_depth=6))
# [3, 7, 15, 1, 292, 1]

# Iterated simultaneous approximation
result = adrius.illl([np.sqrt(2), np.pi], max_denominator=10_000)
for k in range(result.iterations):
    q, p1, p2 = result.relations[k]
    print(f"q={q}  p=({p1},{p2})  quality={result.quality[k]:.2e}")

# LLL lattice reduction
basis = np.array([[1, -1, 3],
                  [1,  0, 5],
                  [1,  2, 6]], dtype=float)
r = adrius.lll_reduce(basis)
print(np.allclose(basis @ r.transform, r.reduced_basis))  # True
```

---

## API reference

### `illl(alpha, *, max_iterations=30, max_denominator=1_000_000, epsilon=0.5, quality_tol=1e-12) → ILLLResult`

Iterated LLL simultaneous Diophantine approximation (Bosma & Smeets 2010).
Produces a sequence of integer relations `(q, p₁, …, pₙ)` satisfying
`max_i |q·αᵢ − pᵢ|` → 0 as `q` grows.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `alpha` | array-like of float | — | Target reals α₁, …, αₙ |
| `max_iterations` | int | 30 | Maximum outer ILLL iterations |
| `max_denominator` | int | 1 000 000 | Stop when q ≥ this value; 0 = no limit |
| `epsilon` | float | 0.5 | Initial lattice scale c₀ = εⁿ⁺¹; must be in (0, 1) |
| `quality_tol` | float | 1e-12 | Stop when max|q·αᵢ − pᵢ| < this value |

Raises `DomainError` if `alpha` is empty or `epsilon` ∉ (0, 1).

---

### `lll_reduce(basis, *, delta=0.75, eta=0.51, max_iter=0) → LLLResult`

LLL lattice basis reduction.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basis` | `np.ndarray`, shape (n, n), float | — | Input basis — **columns** are basis vectors |
| `delta` | float | 0.75 | Lovász condition parameter δ ∈ (0.25, 1.0) |
| `eta` | float | 0.51 | Size-reduction threshold (theoretical min 0.5) |
| `max_iter` | int | 0 | Swap cap; 0 = unlimited |

The transform relationship `basis @ result.transform == result.reduced_basis`
holds exactly up to floating-point noise (‖error‖ < 1e-12 for well-conditioned
inputs).

Raises `DomainError` if the basis is linearly dependent.
Raises `ConvergenceError` if `max_iter` swaps are exceeded.

---

### `cf_expansion(x, *, max_depth=64) → list[int]`

Continued-fraction partial quotients `[a₀; a₁, a₂, …]` of `x`.

---

### `cf_convergents(x, *, max_depth=64) → list[tuple[int, int]]`

Convergents `(pₖ, qₖ)` of `x` for k = 0, 1, 2, …
Each satisfies `|x − pₖ/qₖ| < 1/qₖ²`.

---

### `best_rational(x, *, max_denominator=1_000_000) → Rational`

Best rational approximation `p/q` to `x` with `q ≤ max_denominator`.
Guaranteed optimal by the Hardy-Wright theorem (no `p'/q'` with `q' ≤ max_denominator`
is closer to `x`).

Raises `DomainError` if `max_denominator < 1`.

---

### `simultaneous_approx(alpha, scale, *, delta=0.75, eta=0.51) → SimultApproxResult`

One-shot simultaneous Diophantine approximation via a single LLL call.
Finds integers `q, p₁, …, pₙ` such that `max_i |q·αᵢ − pᵢ| ~ 1/scale^{1/n}`.

> **Note**: For large `scale`, LLL may select an identity column, returning
> `q = 0`.  Use `illl()` for an iterated sequence with guaranteed improving
> quality.

Raises `DomainError` if `alpha` is empty or `scale ≤ 0`.

---

## Result types

### `ILLLResult`

```
.relations   list[list[int]]   [[q, p₁, …, pₙ], …], one per iteration
.quality     list[float]       max_i |q·αᵢ − pᵢ| per iteration
.iterations  int               number of relations found
```

### `LLLResult`

```
.reduced_basis  np.ndarray (float64, n×n)  LLL-reduced basis columns
.transform      np.ndarray (int64,   n×n)  unimodular U: reduced = basis @ U
.swap_count     int                         Lovász swaps performed
```

### `SimultApproxResult`

```
.denominator  int        common denominator q
.numerators   list[int]  [p₁, …, pₙ]
.quality      float      max_i |q·αᵢ − pᵢ|
```

### `Rational`

Exact rational number stored in lowest terms with positive denominator.

```
.numerator    int
.denominator  int
.to_float()   float
.is_integer() bool
.is_zero()    bool
```

Arithmetic operators (`+`, `-`, `*`, `/`, unary `-`) and comparisons (`==`,
`<`, `<=`, `>`, `>=`) are all supported.  `float(r)` and `int(r)` work; `int(r)`
raises `ValueError` if the value is not an integer.

---

## Exception hierarchy

All errors derive from `adrius.Error` which itself derives from
`RuntimeError`.

```python
adrius.Error              # base — also a RuntimeError
├── adrius.DomainError    # bad input (empty alpha, epsilon out of range, etc.)
├── adrius.ConvergenceError  # algorithm hit its iteration cap
└── adrius.PrecisionError    # floating-point precision insufficient
```

```python
try:
    r = adrius.illl([])
except adrius.DomainError as e:
    print(e)          # alpha must be non-empty

# Catch by base class:
except adrius.Error:
    ...
except RuntimeError:
    ...
```

---

## Running the tests

The test suite uses pytest (install with `pip install pytest numpy`).

```bash
# Build the extension first (see Build section above), then:
python -m pytest python/tests/ --rootdir=python -v
```

`python/conftest.py` automatically discovers the compiled extension in any
`build*/python/` directory, so no `PYTHONPATH` setup is needed.

---

## Input conventions

- **`alpha` arrays**: any 1-D array-like — Python lists, NumPy arrays
  (contiguous or not), or other sequences convertible to `float64`.
  Non-contiguous arrays are force-cast to a temporary contiguous buffer.
- **`lll_reduce` basis**: a 2-D NumPy array with `dtype=float`; row-major
  (C-order) arrays are handled correctly by pybind11's Eigen integration.
  **Columns** are the basis vectors (Eigen column-major convention).

---

## References

**Bosma, W., Smeets, I.** (2010). *Finding Good Approximations to Algebraic
Numbers.* — the algorithm underlying `illl()`.

**Lenstra, A.K., Lenstra, H.W., Lovász, L.** (1982). *Factoring polynomials
with rational coefficients.* Mathematische Annalen, 261(4), 515–534. — LLL.

**Hardy, G.H., Wright, E.M.** — *An Introduction to the Theory of Numbers.*
Theorem 171 — `best_rational()` optimality guarantee.
