// Copyright (c) 2025 InGenifold Research LLC. MIT License.
//
// Example: ILLL — Iterated LLL simultaneous Diophantine approximation
//         (Bosma & Smeets, arXiv:1001.4455, 2010)
//
// Usage:  illl <alpha_1> [alpha_2 …] [--max-iter N] [--max-denom D]
//
// Produces a sequence of (q, p₁, …, pₙ) with improving approximation quality
// max_i |q·αᵢ − pᵢ| across iterations.  The Dirichlet coefficient
// q^{1/n} · max|q·αᵢ − pᵢ| is bounded by a constant depending only on n.

#include <adrius/adrius.hpp>

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
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

static std::size_t parse_size(std::string_view s) {
    std::size_t v{};
    auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (ec != std::errc{}) {
        std::cerr << "error: cannot parse '" << s << "' as an integer\n";
        std::exit(1);
    }
    return v;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "usage: illl <alpha_1> [alpha_2 …] [--max-iter N] [--max-denom D]\n";
        return 1;
    }

    std::vector<double> alpha;
    adrius::ILLLParams params{};

    for (int i = 1; i < argc; ++i) {
        std::string_view arg{argv[i]};
        if (arg == "--max-iter") {
            if (i + 1 >= argc) { std::cerr << "error: --max-iter requires a value\n"; return 1; }
            params.max_iterations = parse_size(argv[++i]);
        } else if (arg == "--max-denom") {
            if (i + 1 >= argc) { std::cerr << "error: --max-denom requires a value\n"; return 1; }
            params.max_denominator = static_cast<std::int64_t>(parse_size(argv[++i]));
        } else {
            alpha.push_back(parse_double(arg));
        }
    }

    if (alpha.empty()) {
        std::cerr << "error: at least one alpha value is required\n";
        return 1;
    }

    const std::size_t n = alpha.size();
    std::cout << "Input α: [";
    for (std::size_t i = 0; i < n; ++i)
        std::cout << (i ? ", " : "") << alpha[i];
    std::cout << "]\n\n";

    try {
        auto result = adrius::illl<adrius::EigenBackend>(
            std::span<const double>{alpha}, params);

        // Column widths
        const int w_iter = 5, w_q = 14, w_p = 14, w_qual = 16, w_dc = 18;
        const int total = w_iter + w_q + static_cast<int>(n) * w_p + w_qual + w_dc;

        std::cout << std::setw(w_iter) << "iter"
                  << std::setw(w_q)   << "q";
        for (std::size_t i = 0; i < n; ++i)
            std::cout << std::setw(w_p) << ("p_" + std::to_string(i + 1));
        std::cout << std::setw(w_qual) << "max|q·αᵢ−pᵢ|"
                  << std::setw(w_dc)   << "q^(1/n)·quality"
                  << '\n' << std::string(static_cast<std::size_t>(total), '-') << '\n';

        for (std::size_t k = 0; k < result.relations.size(); ++k) {
            const auto& rel  = result.relations[k];
            const double q   = static_cast<double>(rel[0]);
            const double qual = static_cast<double>(result.quality[k]);
            const double dirichlet = std::pow(q, 1.0 / static_cast<double>(n)) * qual;

            std::cout << std::setw(w_iter) << k + 1
                      << std::setw(w_q)    << rel[0];
            for (std::size_t i = 0; i < n; ++i)
                std::cout << std::setw(w_p) << rel[i + 1];
            std::cout << std::setw(w_qual) << std::scientific << std::setprecision(4) << qual
                      << std::setw(w_dc)   << std::fixed     << std::setprecision(6) << dirichlet
                      << '\n';
        }

        std::cout << '\n' << result.iterations << " approximation(s) found.\n";
    } catch (const adrius::Error& e) {
        std::cerr << "adrius error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
