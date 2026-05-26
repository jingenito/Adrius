// Copyright (c) 2025 InGenifold Research LLC. MIT License.
//
// Example: LLL lattice basis reduction
//
// Reduces a "bad" 3×3 basis using the two-stage API and prints the
// reduced basis, unimodular transform U (reduced = original · U),
// and the number of Lovász swaps performed.

#include <adrius/adrius.hpp>

#include <iomanip>
#include <iostream>

int main() {
    // Classic example from LLL (1982): three nearly-parallel vectors.
    adrius::EigenBackend::matrix_type basis(3, 3);
    basis << 1, -1,  3,
             1,  0,  5,
             1,  2,  6;

    std::cout << "Input basis (columns):\n" << basis << "\n\n";

    // Stage 1: validate and compute initial Gram-Schmidt decomposition.
    auto prepared = adrius::preprocess_lll<adrius::EigenBackend>(basis);

    // Stage 2: iterate with Lovász incremental updates.
    auto result = adrius::lll_reduce<adrius::EigenBackend>(std::move(prepared));

    std::cout << "Reduced basis:\n" << result.reduced_basis << "\n\n";
    std::cout << "Unimodular transform U  (reduced = original · U):\n"
              << result.transform << "\n\n";
    std::cout << "Lovász swaps: " << result.swap_count << "\n";

    // Verify the transform relationship B' = B · U.
    const double err =
        (basis * result.transform.cast<double>() - result.reduced_basis).norm();
    std::cout << "Verification ||B·U − B'|| = " << std::scientific << err << "\n";

    return 0;
}
