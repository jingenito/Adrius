# Contributing to Adrius

Thank you for your interest in contributing to Adrius! This document explains how to get started.

## Types of Contributions Welcome

- **New algorithms**: Simultaneous Diophantine approximation, lattice reduction variants, etc.
- **Algorithm improvements**: Better bounds, faster convergence, numerical stability
- **Backends**: New linear algebra backends (e.g., Armadillo, Blaze) beyond Eigen
- **Tests**: Property-based tests, edge case coverage, benchmarks
- **Documentation**: Algorithm explanations, usage examples, design rationale
- **Bug fixes**: Correctness issues, numerical problems, platform-specific issues

**Note**: Major architectural changes or new data structures should be discussed as an issue first.

## Getting Started

### 1. Read the Project Documentation
- **CLAUDE.md** — Code style rules and design principles
- **DEVELOPMENT.md** — Development setup, build instructions, and common workflows
- **docs/design-decisions.md** — Architectural decisions (why the library is structured this way)
- **README.md** — User-facing overview

### 2. Set Up Your Environment
Follow the "Quick Start" section of **DEVELOPMENT.md**:
```bash
git clone https://github.com/jingenito/Adrius.git
cd Adrius
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

### 3. Create a Feature Branch
Always branch off the latest `main`:
```bash
git checkout main
git pull origin main
git checkout -b feature/your-feature-name
```

See "Branch Naming" in CLAUDE.md for naming conventions.

## Development Workflow

### Before You Start Coding
1. **Algorithm changes**: Read the paper and relevant algorithm docs first
2. **New algorithm**: Check DEVELOPMENT.md "Adding a New Algorithm" section
3. **Bug fixes**: Understand the root cause before implementing a fix

### While Coding
1. **Follow CLAUDE.md** style rules (code style, API design, memory model, etc.)
2. **Write tests first** (or alongside code): At minimum, hand-computed examples
3. **Document as you go**: Comments explain "why", not "what"
4. **Run tests frequently**: `ctest --test-dir build --output-on-failure`

### Before Creating a Pull Request
1. **Run full test suite**: `ctest --test-dir build --output-on-failure`
2. **All tests must pass**: Fix failures before submitting PR
3. **No compiler warnings**: Check build output for W4 warnings (MSVC) or -Wall (GCC/Clang)
4. **Document algorithm changes**: 
   - Create `docs/algorithm-name.md` for new algorithms
   - Include research paper references
   - Add examples to README.md
5. **Update CLAUDE.md or DEVELOPMENT.md** only if adding new guidelines

## Pull Request Process

### Creating a PR
1. Push your branch: `git push -u origin feature/your-feature-name`
2. Create PR on GitHub
3. Use the PR template (`.github/pull_request_template.md`):
   - Summary of changes
   - Type of change (Feature/Fix/Docs/etc.)
   - Test results (copy output from `ctest`)
   - **Research papers referenced** (with specific sections/equations)
   - Checklist items completed

### PR Description Example
```markdown
## Summary
Implement ILLL algorithm from Bosma & Smeets (2010) for finding simultaneous 
Diophantine approximations with prescribed quality bounds.

## Type of Change
- [x] Feature

## Research Papers
- **Bosma, W. & Smeets, I.** (2010). "Finding simultaneous Diophantine 
  approximations with prescribed quality." arXiv:1001.4455.
  - Sections: Algorithm 1, Lemma 3.4, Theorem 3.5 (eq. 24, 28)

## Testing
All tests pass locally:
```
100% tests passed, 0 tests failed out of 46
Total Test time (real) = 0.53 sec
```

## Documentation
- [x] New algorithm documentation in `docs/illl-algorithm.md`
- [x] README.md updated with ILLL example
- [x] Test suite covers hand-computed examples and bounds
```

### Code Review
- Reviews will check:
  - Correctness against referenced papers
  - Code style compliance with CLAUDE.md
  - Test coverage
  - Documentation completeness
  - Performance and numerical stability
- Address feedback with new commits (don't rebase/force-push unless requested)

### Merging
- Once approved, PR will be merged to `main`
- Feature branch is deleted automatically
- Commit history is preserved (not squashed)

## Code Style Guide

### Must-Follow Rules
- **No `using namespace` at file scope** (use `using` in function scope or qualify names)
- **`[[nodiscard]]`** on every function returning a result struct or prepared state
- **Named-field structs** for multi-value returns (never `std::tuple`)
- **Two-stage API**: `preprocess_foo()` → `PreparedFoo`, then `foo(prepared)`
- **Move-only state**: `PreparedFoo` types disable copy (only move)
- **Copyright header**: First line of every file is `// Copyright (c) 2025 InGenifold Research LLC. MIT License.`

### Commit Message Format
```
Implement new algorithm from Paper et al. (YYYY)

- Detailed bullet points of what was added
- Explain key design decisions
- Reference paper sections if relevant

Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>
```

### Comments
- **Good**: "Rescale only column 0; keeps basis nearly reduced for amortized cost"
- **Bad**: "Scale column 0"
- **Good**: "Lovász incremental update (LLL paper §3) avoids full recompute after swap"
- **Bad**: "Update Gram-Schmidt vectors"

## Algorithm Implementation Checklist

For new algorithms, ensure:
- [ ] Two-stage API (preprocess + iterate)
- [ ] Result struct with `[[nodiscard]]` and named fields
- [ ] Prepared state struct (move-only)
- [ ] Params struct with sensible defaults
- [ ] At least one hand-computed example test
- [ ] Tests for correctness and mathematical bounds
- [ ] Algorithm documentation in `docs/algorithm-name.md`
- [ ] Usage examples in README.md
- [ ] Research papers cited in docs and PR description
- [ ] All tests passing
- [ ] No compiler warnings

## Testing Guidelines

### Test Names
```cpp
// Good: describes what's verified
TEST(ILLL, OneDimensionalBoundedDirichletCoefficient)
TEST(Rational, NegativeDenominatorNormalized)

// Bad: generic names
TEST(ILLL, Test1)
TEST(Rational, BasicTest)
```

### Hand-Computed Examples
Every algorithm must have at least one test with values you can verify by hand:
```cpp
TEST(MyAlgorithm, HandComputedExample) {
    const auto result = adrius::my_algorithm<B>({1.41, 3.14});
    EXPECT_EQ(result.denominator, 70);  // Known from manual calculation
}
```

### Bounds Testing
Verify algorithms satisfy their mathematical guarantees:
```cpp
TEST(ILLL, BosmaSmeetsBound) {
    const auto result = adrius::illl<B>(alpha, params);
    for (const auto& rel : result.relations) {
        const double quality = /* compute max approximation error */;
        EXPECT_LT(rel.denominator * quality, UPPER_BOUND);
    }
}
```

## Questions or Issues?

- **Code style questions**: Check CLAUDE.md
- **Setup/build issues**: Check DEVELOPMENT.md "Troubleshooting"
- **Algorithm questions**: Check relevant `docs/algorithm-name.md`
- **Unsure about approach**: Open an issue to discuss before coding

## License

All contributions are made under the MIT License. By contributing, you agree that your contributions may be used under this license.

---

**Happy coding! 🚀**
