# Adrius Development Guide

This guide covers environment setup, building, testing, and common development workflows for Adrius contributors.

## Prerequisites

### Compiler & Tools
- **C++20 compiler**: MSVC 19.29+, GCC 11+, or Clang 13+
- **CMake**: 3.21 or later
- **Git**: For version control

### Platform-Specific Setup

#### Windows (MSVC)
```bash
# Install via winget (requires Windows 11 or winget on Windows 10)
winget install cmake
winget install --id GitHub.cli
```

If using Visual Studio Build Tools manually:
- Download from https://visualstudio.microsoft.com/downloads/ (select "Build Tools for Visual Studio")
- Install C++ workload
- CMake will auto-detect cl.exe compiler

#### Linux (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install build-essential cmake git
```

#### macOS
```bash
# Install Xcode Command Line Tools
xcode-select --install

# Install CMake via Homebrew
brew install cmake
```

## Quick Start

### 1. Clone and Setup
```bash
git clone https://github.com/jingenito/Adrius.git
cd Adrius
```

### 2. Configure Build
```bash
# Linux / macOS
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

# Windows (with VS 2022 Build Tools)
# Note: On Windows, you may need to run vcvarsall.bat first to set up MSVC environment
cmake -B build --preset windows-msvc-debug
```

CMake will automatically fetch dependencies (Eigen3, GoogleTest).

### 3. Build
```bash
cmake --build build

# Or specify configuration on Windows
cmake --build build --config Debug
```

### 4. Run Tests
```bash
ctest --test-dir build --output-on-failure
```

## Development Workflows

### Adding a New Algorithm

1. **Create feature branch** off latest main:
   ```bash
   git checkout main
   git pull origin main
   git checkout -b feature/my-algorithm-name
   ```

2. **Create header in appropriate directory**:
   - `include/adrius/approx/my_algorithm.hpp` for approximation algorithms
   - `include/adrius/linalg/my_algorithm.hpp` for linear algebra
   - `include/adrius/util/my_type.hpp` for utility types

3. **Implement two-stage API**:
   ```cpp
   template <Backend B = DefaultBackend>
   [[nodiscard]] struct MyAlgorithmResult {
       // Named fields for results
       std::vector<result_type> results;
       std::size_t iterations{0};
   };

   template <Backend B = DefaultBackend>
   struct PreparedMyAlgorithm {
       // Move-only state
       matrix_of<B> prepared_matrix;
       // No default constructor, only move
       PreparedMyAlgorithm() = delete;
       PreparedMyAlgorithm(PreparedMyAlgorithm&&) noexcept = default;
       PreparedMyAlgorithm& operator=(PreparedMyAlgorithm&&) noexcept = default;
       PreparedMyAlgorithm(const PreparedMyAlgorithm&) = delete;
       PreparedMyAlgorithm& operator=(const PreparedMyAlgorithm&) = delete;
   };

   // Stage 1: Preprocessing
   template <Backend B = DefaultBackend>
   [[nodiscard]] PreparedMyAlgorithm<B>
   preprocess_my_algorithm(/* inputs */, MyAlgorithmParams params = {}) { /* ... */ }

   // Stage 2: Iteration
   template <Backend B = DefaultBackend>
   [[nodiscard]] MyAlgorithmResult<B>
   my_algorithm(PreparedMyAlgorithm<B> prepared, MyAlgorithmParams params = {}) { /* ... */ }

   // Convenience overload: combines both stages
   template <Backend B = DefaultBackend>
   [[nodiscard]] MyAlgorithmResult<B>
   my_algorithm(/* raw inputs */, MyAlgorithmParams params = {}) {
       return my_algorithm<B>(preprocess_my_algorithm<B>(/* inputs */, params), params);
   }
   ```

4. **Add params struct** to `include/adrius/core/params.hpp`:
   ```cpp
   struct MyAlgorithmParams {
       std::size_t max_iterations = 30;
       double tolerance = 1e-12;
       // ... with sensible defaults
   };
   ```

5. **Write comprehensive tests** in `tests/test_my_algorithm.cpp`:
   - At least one hand-computed small example
   - Tests for correctness and bounds
   - Edge cases and error conditions
   - Use descriptive test names

6. **Create algorithm documentation** in `docs/my_algorithm.md`:
   - Overview and mathematical foundation
   - Implementation details with paper references
   - Usage examples
   - Design decisions
   - References section

7. **Update README.md** if adding public API

8. **Commit and push**:
   ```bash
   git add include/adrius/approx/my_algorithm.hpp tests/test_my_algorithm.cpp docs/my_algorithm.md
   git commit -m "Implement my algorithm from Paper et al. (YYYY)

   - Add MyAlgorithmResult<B> struct
   - Implement preprocess_my_algorithm<B>() and my_algorithm<B>()
   - Add comprehensive test suite (N tests)
   - Add algorithm documentation with paper references

   Co-Authored-By: Claude Haiku 4.5 <noreply@anthropic.com>
   "
   git push -u origin feature/my-algorithm-name
   ```

### Running Specific Tests

```bash
# All tests
ctest --test-dir build --output-on-failure

# Only ILLL tests
ctest --test-dir build --output-on-failure -R ILLL

# Only Rational tests
ctest --test-dir build --output-on-failure -R Rational

# One specific test
ctest --test-dir build --output-on-failure -R "Rational.Addition"

# Run test executable directly with gtest filter
./build/tests/test_illl --gtest_filter="ILLL.*"
```

### Rebuilding After Header Changes

Header-only code requires rebuilds when includes change. If cmake isn't picking up changes:

```bash
# Touch test file to force rebuild
touch tests/test_name.cpp

# Rebuild specific target
cmake --build build --target test_name
```

### Debugging Tests

```bash
# Run test with verbose output
ctest --test-dir build --output-on-failure -VV -R "TestName"

# Run test executable directly (gives full gtest output)
./build/tests/test_name.exe --gtest_filter="Suite.TestName"

# Run with gdb/lldb (Linux/macOS)
gdb ./build/tests/test_name
```

### Checking Code Style

While Adrius doesn't enforce automatic formatting, follow the rules in CLAUDE.md:
- No `using namespace` at file scope
- `[[nodiscard]]` on result structs and prepared state
- Named structs instead of tuples
- Copyright header on every file
- Comments explain "why", not "what"

Before submitting PR, do a manual review:
```bash
git diff main...HEAD  # Show all changes in this branch
```

## Build Configurations

### Development (Debug)
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

### Release
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

### With System Eigen (skip auto-fetch)
```bash
cmake -B build -DADRIUS_USE_SYSTEM_EIGEN=ON
```

### With Boost.Multiprecision Backend
```bash
# Auto-fetch Boost
cmake -B build -DADRIUS_ENABLE_BOOST_MULTIPRECISION=ON

# Or use system Boost (if installed)
cmake -B build -DADRIUS_ENABLE_BOOST_MULTIPRECISION=ON -DADRIUS_USE_SYSTEM_BOOST=ON
```

Note: First build may take 2-5 minutes to build Boost dependencies. See `docs/boost-multiprecision-backend.md` for usage.

### Disable Tests
```bash
cmake -B build -DADRIUS_BUILD_TESTS=OFF
```

### Disable Examples
```bash
cmake -B build -DADRIUS_BUILD_EXAMPLES=OFF
```

## Troubleshooting

### CMake Can't Find Compiler (Windows MSVC)

**Problem**: CMake says "No CMAKE_CXX_COMPILER could be found"

**Solution**: Ensure MSVC environment is initialized:
```bash
# In PowerShell (adjust path for your VS version)
$vsPath = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
& "$vsPath\VC\Auxiliary\Build\vcvarsall.bat" x64
# Then run cmake
```

Or use CMake GUI to select compiler explicitly.

### IntelliSense Errors but Code Compiles

**Problem**: VS Code reports missing headers like `<concepts>` even though code builds

**Solution**: 
1. Ensure `.vscode/c_cpp_properties.json` has `compileCommands` pointing to `build/compile_commands.json`
2. Add CMake configure flag: `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`
3. Reload VS Code C++ extension

### Test Executable Won't Run (Windows)

**Problem**: Linker errors or "cannot find -lgtest" on Windows

**Solution**: Ensure GoogleTest was built:
```bash
cmake --build build  # Builds all targets including gtest
ctest --test-dir build --output-on-failure
```

### Eigen Compilation Warnings

**Problem**: MSVC C5054 deprecation warnings from Eigen headers

**Solution**: Already handled by pragma guards in `backend/eigen.hpp`. If you see warnings, check that file includes the pragma suppression.

## Continuous Integration

All commits pushed to `main` run CI (GitHub Actions). Before pushing, ensure:
1. All tests pass: `ctest --test-dir build --output-on-failure`
2. No compiler errors or warnings
3. Code follows style rules in CLAUDE.md

If CI fails on your PR:
1. Check the error message in GitHub Actions logs
2. Fix locally: `git add . && git commit ...`
3. Push fix: `git push`
4. CI will re-run automatically

## Documentation

- **CLAUDE.md** — Code style and design rules (check here first!)
- **docs/design-decisions.md** — Architectural choices for the library
- **docs/illl-algorithm.md** — ILLL algorithm documentation (reference for implementing similar algorithms)
- **README.md** — User-facing overview and usage examples
- **DEVELOPMENT.md** — This file

For algorithm documentation, model your `docs/algorithm-name.md` after `docs/illl-algorithm.md`.

## Getting Help

- Check CLAUDE.md for code style questions
- Check docs/design-decisions.md for architectural questions
- Check individual algorithm docs (e.g., docs/illl-algorithm.md) for implementation questions
- Run specific tests in verbose mode to debug: `ctest --test-dir build -VV -R TestName`
