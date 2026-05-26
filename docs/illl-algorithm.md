# ILLL: Iterated LLL for Simultaneous Diophantine Approximation

This document describes the Iterated LLL (ILLL) algorithm as implemented in Adrius, including the mathematical foundation, implementation details, and usage patterns.

## Overview

ILLL is an iterative algorithm that finds simultaneous Diophantine approximations to n real numbers with bounded quality. Given targets α = (α₁, α₂, ..., αₙ), the algorithm produces a sequence of solutions:

- (q₁, p₁,₁, p₁,₂, ..., p₁,ₙ)
- (q₂, p₂,₁, p₂,₂, ..., p₂,ₙ)
- ...

where each solution satisfies: **|q·αᵢ − pᵢ| ≤ quality** for all i.

The denominators q grow monotonically, and the quality (approximation error) improves geometrically per iteration.

## Algorithm Components

### 1. Lattice Preprocessing

The algorithm begins by constructing an (n+1)×(n+1) lattice matrix:

```
L₀ = [c₀   0   0  ...  0  ]
     [α₁   1   0  ...  0  ]
     [α₂   0   1  ...  0  ]
     [...  ... ... ... ... ]
     [αₙ   0   0  ...  1  ]
```

where:
- **c₀ = ε^{n+1}** with ε ∈ (0,1), default ε = 0.5
- For n = 2: c₀ = 0.125
- **Column 0** encodes the denominator (scaled by c₀) and approximation targets
- **Columns 1..n** form the identity block

#### Why This Construction?

A naive lattice (N, Nα₁, ..., Nαₙ) fails because:
- First column has norm ≈ N, which is much larger than the identity columns (norm ≈ 1)
- LLL's reduction algorithm minimizes norms, so it ignores the first column entirely
- Result: q → 0 (no useful denominator found)

The Bosma-Smeets construction balances all column norms at O(1), forcing LLL to treat the denominator direction as equally important as the error directions.

### 2. Main Iteration

For k = 0, 1, 2, ..., the algorithm:

#### Step 1: LLL Reduction
Apply LLL reduction to the current basis Lₖ, producing a reduced basis B.

#### Step 2: Extract Approximation
From the first column of B:
- **b₀** = first element of first column
- **bᵢ** = (i+1)-th element, for i = 1..n

Compute:
```
q = round(|b₀| / cₖ)           // denominator
pᵢ = round(q·αᵢ − bᵢ)           // numerators
quality = max_i |bᵢ|           // approximation error
```

If q > 0:
- Record the relation (q, p₁, ..., pₙ)
- Check termination conditions (see below)

#### Step 3: Scale Update
Update the scale factor:
```
cₖ₊₁ = cₖ · 2^{−n(n+1)/4}
```

This comes from Bosma & Smeets Lemma 3.4: the scale decays geometrically such that per-iteration quality improves by a constant factor depending only on n.

#### Step 4: Rescale Basis (warm start)
`c_k` appears exclusively in row 0 of the lattice matrix.  Transitioning
from Λ_k to Λ_{k+1} therefore requires scaling **row 0 of every column** by
`c_{k+1}/c_k`:
```
Lₖ₊₁[0, j] ← B[0, j] · (cₖ₊₁ / cₖ)   for all j = 0..n
Lₖ₊₁[i, j] ← B[i, j]                  for i ≥ 1
```

Rows 1..n encode ℤ-linear combinations of the αᵢ, which are unchanged when
the scale c_k changes.  Scaling column 0 (a common error) corrupts those rows
and locks the extracted denominator at q = 1 for all subsequent iterations.

This warm start leaves the basis nearly reduced for Λ_{k+1}, amortizing the
cost of repeated LLL calls.

### 3. Termination

The algorithm stops when ANY condition is met:

1. **Denominator cap**: q ≥ max_denominator
2. **Quality threshold**: max_i |qαᵢ − pᵢ| < quality_tol
3. **Numerical underflow**: cₖ < 1e-14 (prevents gram_schmidt linear dependence)
4. **Iteration limit**: k ≥ max_iterations

## Mathematical Guarantees

### Theorem (Bosma & Smeets 3.5)
The Dirichlet coefficient satisfies:
```
θ(k) := q(k)^{1/n} · max_i |q(k)·αᵢ − p(k,i)| ≤ C
```
where C is a constant depending only on n (not on the specific targets α).

### Convergence
The algorithm finds all approximations up to denominator q_max in **O(log q_max)** iterations.

Per-iteration quality decay:
```
quality(k) ≈ quality(0) · 2^{−kn/4}
```

## Implementation Details

### Scale Formula Derivation

From the paper (eq. 28), with d = 2 (reduction factor in log₂ base):

```
cₖ = (2^{−kn/4} · ε)^{n+1}

Scale ratio:
cₖ₊₁ / cₖ = (2^{−(k+1)n/4} / 2^{−kn/4})^{n+1}
           = 2^{−n/4·(n+1)}
           = 2^{−n(n+1)/4}
```

### Numerical Safeguards

**Guard condition (added, not in paper):**
```cpp
if (iter > 0 && prepared.scale * prepared.scale < params.lll.gso.zero_threshold)
    break;
```

Comparing `c_k²` against `zero_threshold` (a squared-norm threshold) is dimensionally consistent.  When `c_k²` falls below the Gram-Schmidt linear-dependence threshold, column 0's squared norm would be indistinguishable from zero and LLL would throw `DomainError`.

### Two-Stage API

The algorithm is split into two phases for flexibility:

```cpp
// Phase 1: Preprocessing
auto prepared = adrius::preprocess_illl<B>(alpha, params);

// Phase 2: Iteration
auto result = adrius::illl<B>(prepared, params);
```

Benefits:
- Callers can inspect or cache the prepared lattice
- Multiple invocations can share preprocessing (e.g., parameter sweeps)
- Clear separation of concerns

## Usage Example

```cpp
#include <adrius/adrius.hpp>
#include <iostream>
#include <numbers>
#include <span>

int main() {
    // Find simultaneous approximations to √2 and π
    const std::vector<double> alpha = {std::numbers::sqrt2, std::numbers::pi};
    
    adrius::ILLLParams params;
    params.max_iterations  = 20;
    params.max_denominator = 100'000;
    params.quality_tol     = 1e-10;
    params.epsilon         = 0.5;  // c₀ = 0.5³ = 0.125 for n=2
    
    auto result = adrius::illl<adrius::EigenBackend>(
        std::span<const double>{alpha}, 
        params
    );
    
    for (std::size_t k = 0; k < result.iterations; ++k) {
        const auto& rel = result.relations[k];
        std::cout << "q = " << rel[0] 
                  << ", p₁ = " << rel[1] 
                  << ", p₂ = " << rel[2]
                  << ", quality = " << result.quality[k] << '\n';
    }
    
    return 0;
}
```

## Key Design Decisions

1. **Balanced Lattice**: c₀ = ε^{n+1} ensures all columns have comparable norm, not ignoring denominator direction.

2. **Geometric Scale Decay**: 2^{−n(n+1)/4} per iteration ensures quality improves predictably and Bosma-Smeets bounds hold.

3. **Row-0 Rescaling**: Row 0 of all columns is rescaled between iterations (not just column 0). This correctly maps Λ_k → Λ_{k+1} while preserving the ℤ-linear error rows, keeping the basis nearly reduced and amortizing LLL cost.

4. **Move-Only State**: PreparedILLL<B> is move-only to prevent accidental state sharing and ensure proper cleanup.

5. **Two-Stage API**: Preprocessing and iteration are separate, following Adrius design principles for composability and introspection.

6. **Numerical Safeguards**: Check `c_k² < zero_threshold` (squared, to match the GSO's squared-norm comparisons) before LLL to prevent floating-point underflow errors.

## References

**Bosma, W. & Smeets, I.** (2010). *Finding simultaneous Diophantine approximations with prescribed quality.* arXiv:1001.4455.

Key sections:
- §2.1: Lattice construction and initialization
- Lemma 3.4: Scale decay formula
- Theorem 3.5: Dirichlet coefficient bound and convergence guarantee
- Algorithm 1: Main iteration framework
- Equation (28): Scale update formula

The implementation follows the paper precisely, with only the addition of a numerical safeguard (scale < 1e-14) for double-precision stability.

## Testing

The algorithm is tested against:
- **Correctness**: Known convergents of √2 are recovered
- **Bounds**: 1D approximations satisfy Bosma-Smeets Dirichlet coefficient bounds
- **Monotonicity**: Denominators increase monotonically across iterations
- **Termination**: Algorithm respects max_denominator and quality_tol caps
- **Numerical Stability**: No crash on ill-conditioned inputs

See `tests/test_illl.cpp` for comprehensive test suite (8 tests covering all aspects).
