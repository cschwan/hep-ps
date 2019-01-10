#include "hep/ps/phase_space_tools.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <utility>

namespace hep
{

template <typename T>
void print_summary_of_invariants(std::vector<T> const& phase_space)
{
    using std::fabs;

    std::cout << "PHASE SPACE POINT\n=================\n";

    std::size_t const n = phase_space.size() / 4;

    for (std::size_t i = 0; i != n; ++i)
    {
        std::cout << 'P' << i << " = {"
            << phase_space.at(4 * i + 0) << ", " << phase_space.at(4 * i + 1) << ", "
            << phase_space.at(4 * i + 2) << ", " << phase_space.at(4 * i + 3) << "}\n";
    }

    std::cout << '\n';

    std::vector<std::pair<std::size_t, T>> invariants;
    invariants.reserve(1 << n);

    for (std::size_t i = 1; i != invariants.capacity() + 1; ++i)
    {
        std::array<T, 4> sum = {};

        for (std::size_t j = 0; j != n; ++j)
        {
            if ((i & (1 << j)) != 0)
            {
                T const sign = (j < 2) ? T(-1.0) : T(1.0);

                sum[0] += sign * phase_space[4 * j + 0];
                sum[1] += sign * phase_space[4 * j + 1];
                sum[2] += sign * phase_space[4 * j + 2];
                sum[3] += sign * phase_space[4 * j + 3];
            }
        }

        invariants.emplace_back(i, sum[0]*sum[0] - sum[1]*sum[1] - sum[2]*sum[2] - sum[3]*sum[3]);
    }

    std::sort(invariants.begin(), invariants.end(), [](std::pair<std::size_t, T> const& a,
        std::pair<std::size_t, T> const& b) { return fabs(a.second) < fabs(b.second); });

    std::cout << "SUMMARY OF INVARIANTS\n=====================\n"
        << 's' << invariants.at(0).first
        << '=' << invariants.at(0).second << '\n'
        << 's' << invariants.at(1).first
        << '=' << invariants.at(1).second << '\n'
        << 's' << invariants.at(2).first
        << '=' << invariants.at(2).second << '\n'
        << "   ...\n"
        << 's' << invariants.back().first
        << '=' << invariants.back().second << "\n\n";
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template void print_summary_of_invariants(std::vector<double> const&);

}
