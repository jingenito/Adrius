// Copyright (c) 2025 InGenifold Research LLC. MIT License.
//
// Example: best rational approximation via continued fractions
//
// Usage:  best_rational <real> [max_denominator]
//
// Prints the CF convergent table for the given real number and then
// shows the best rational approximation p/q with denominator ≤ max_denominator.

#include <adrius/adrius.hpp>

#include <charconv>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string_view>

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
        std::cerr << "usage: best_rational <real> [max_denominator]\n";
        return 1;
    }

    const double x = parse_double(argv[1]);
    const std::int64_t max_denom = argc >= 3
        ? static_cast<std::int64_t>(std::atoll(argv[2]))
        : std::int64_t{1000};

    // Convergent table
    auto convs = adrius::cf_convergents(x, adrius::CFParams{.max_depth = 20});
    std::cout << "Convergents for " << x << ":\n";
    std::cout << std::setw(4)  << "k"
              << std::setw(14) << "p_k"
              << std::setw(14) << "q_k"
              << std::setw(20) << "p_k / q_k"
              << std::setw(16) << "|x - p/q|"
              << '\n' << std::string(68, '-') << '\n';

    for (std::size_t k = 0; k < convs.size(); ++k) {
        const double approx = static_cast<double>(convs[k].p)
                            / static_cast<double>(convs[k].q);
        std::cout << std::setw(4)  << k
                  << std::setw(14) << convs[k].p
                  << std::setw(14) << convs[k].q
                  << std::setw(20) << std::fixed << std::setprecision(12) << approx
                  << std::setw(16) << std::scientific << std::setprecision(4)
                  << std::abs(x - approx) << '\n';
    }

    // Best rational within the denominator bound
    auto best = adrius::best_rational(x, {.max_denominator = max_denom});
    const double best_approx = best.to_double();
    std::cout << "\nBest approximation with q ≤ " << max_denom << ":  "
              << best
              << "  ≈ " << std::fixed << std::setprecision(12) << best_approx
              << "  (error " << std::scientific << std::setprecision(4)
              << std::abs(x - best_approx) << ")\n";

    return 0;
}
