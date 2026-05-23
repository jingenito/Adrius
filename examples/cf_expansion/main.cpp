// Copyright (c) 2025 InGenifold Research LLC. MIT License.
//
// Example: continued fraction expansion
//
// Usage:  cf_expansion <real> [max_depth]
//
// Prints the partial quotients [a_0; a_1, a_2, …] and the convergents p_k/q_k
// for the given real number.

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
        std::cerr << "usage: cf_expansion <real> [max_depth]\n";
        return 1;
    }

    const double x         = parse_double(argv[1]);
    const std::size_t depth = argc >= 3
        ? static_cast<std::size_t>(std::atoi(argv[2]))
        : std::size_t{32};

    adrius::CFParams params{.max_depth = depth};

    // Print the continued fraction notation
    std::cout << "CF(" << x << ") = [";
    bool first = true;
    for (std::int64_t a : adrius::cf_view(x, params)) {
        if (!first) std::cout << "; ";
        std::cout << a;
        first = false;
    }
    std::cout << "]\n\n";

    // Print convergent table
    auto convs = adrius::cf_convergents(x, params);
    std::cout << std::setw(4) << "k"
              << std::setw(14) << "p_k"
              << std::setw(14) << "q_k"
              << std::setw(20) << "p_k/q_k"
              << std::setw(16) << "|x - p/q|"
              << '\n'
              << std::string(68, '-') << '\n';

    for (std::size_t k = 0; k < convs.size(); ++k) {
        const double approx = static_cast<double>(convs[k].p)
                              / static_cast<double>(convs[k].q);
        std::cout << std::setw(4) << k
                  << std::setw(14) << convs[k].p
                  << std::setw(14) << convs[k].q
                  << std::setw(20) << std::setprecision(12) << approx
                  << std::setw(16) << std::setprecision(4) << std::scientific
                  << std::abs(x - approx)
                  << '\n';
    }

    return 0;
}
