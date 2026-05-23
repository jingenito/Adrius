// Copyright (c) 2025 InGenifold Research LLC. MIT License.
//
// Example: simultaneous Diophantine approximation (Bosma-Smeets workflow)
//
// Usage:  bosma_smeets <alpha_1> <alpha_2> … [--scale N]
//
// Finds a common denominator q and integers p_1, …, p_n such that
// |q·αᵢ − pᵢ| is small for all i simultaneously, using LLL lattice reduction.

#include <adrius/adrius.hpp>

#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <string_view>
#include <vector>

static double parse_double(std::string_view s) {
    double v{};
    auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{}) {
        std::cerr << "error: cannot parse '" << s << "' as a real number\n";
        std::exit(1);
    }
    return v;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "usage: bosma_smeets <alpha_1> [alpha_2 …] [--scale N]\n";
        return 1;
    }

    std::vector<double> alpha;
    double scale = 1e6;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg{argv[i]};
        if (arg == "--scale") {
            if (i + 1 >= argc) { std::cerr << "error: --scale requires a value\n"; return 1; }
            scale = parse_double(argv[++i]);
        } else {
            alpha.push_back(parse_double(arg));
        }
    }

    if (alpha.empty()) {
        std::cerr << "error: at least one alpha value is required\n";
        return 1;
    }

    std::cout << "Input α: [";
    for (std::size_t i = 0; i < alpha.size(); ++i)
        std::cout << (i ? ", " : "") << alpha[i];
    std::cout << "]\nScale N = " << scale << "\n\n";

    try {
        auto result = adrius::simultaneous_approx<adrius::EigenBackend>(
            std::span<const double>{alpha}, scale);

        std::cout << "Denominator q = " << result.denominator << "\n";
        std::cout << "Numerators:    [";
        for (std::size_t i = 0; i < result.numerators.size(); ++i)
            std::cout << (i ? ", " : "") << result.numerators[i];
        std::cout << "]\n";
        std::cout << "Quality (max |q·αᵢ − pᵢ|) = "
                  << std::scientific << std::setprecision(4) << result.quality << "\n\n";

        std::cout << "Approximation detail:\n";
        std::cout << std::setw(6)  << "i"
                  << std::setw(16) << "αᵢ"
                  << std::setw(14) << "pᵢ"
                  << std::setw(20) << "pᵢ/q"
                  << std::setw(16) << "|q·αᵢ − pᵢ|"
                  << '\n' << std::string(72, '-') << '\n';

        for (std::size_t i = 0; i < alpha.size(); ++i) {
            const double approx = static_cast<double>(result.numerators[i])
                                  / static_cast<double>(result.denominator);
            const double err    = std::abs(static_cast<double>(result.denominator) * alpha[i]
                                           - static_cast<double>(result.numerators[i]));
            std::cout << std::setw(6)  << i
                      << std::setw(16) << std::fixed << std::setprecision(10) << alpha[i]
                      << std::setw(14) << result.numerators[i]
                      << std::setw(20) << std::setprecision(10) << approx
                      << std::setw(16) << std::scientific << std::setprecision(4) << err
                      << '\n';
        }
    } catch (const adrius::Error& e) {
        std::cerr << "adrius error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
