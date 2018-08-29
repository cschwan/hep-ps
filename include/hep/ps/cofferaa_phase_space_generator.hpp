#ifndef HEP_PS_COFFERAA_PHASE_GENERATOR_HPP
#define HEP_PS_COFFERAA_PHASE_GENERATOR_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017  Christopher Schwan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "hep/ps/phase_space_generator.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace hep
{

/// A modified version of the phase space generator used in \f$
/// \mathrm{Coffer}\gamma\gamma \f$, see \cite Bredenstein:2005zk .
template <typename T>
std::unique_ptr<phase_space_generator<T>> make_cofferaa_phase_space_generator(
    T min_energy,
    T cmf_energy,
    std::vector<int> const& process,
    lusifer_constants<T> const& constants,
    std::vector<std::tuple<int, int, int>> const& dipoles
        = std::vector<std::tuple<int, int, int>>{},
    std::size_t extra_random_numbers = 0
);

template <typename T>
std::unique_ptr<phase_space_generator<T>>
make_minimal_cofferaa_phase_space_generator(
    T min_energy,
    T cmf_energy,
    std::vector<int> const& process,
    lusifer_constants<T> const& constants,
    std::vector<std::tuple<int, int, int>> const& dipoles
        = std::vector<std::tuple<int, int, int>>{},
    std::size_t extra_random_numbers = 0
);

}

#endif
