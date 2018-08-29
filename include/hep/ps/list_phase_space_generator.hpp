#ifndef HEP_PS_LIST_PHASE_SPACE_GENERATOR_HPP
#define HEP_PS_LIST_PHASE_SPACE_GENERATOR_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2018  Christopher Schwan
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
#include "hep/ps/phase_space_generator.hpp"

#include <memory>
#include <vector>

namespace hep
{

/// Class capturing a phase space point together with the luminosity
/// information.
template <typename T>
class list_phase_space_point
{
public:
    /// Constructor.
    list_phase_space_point(
        luminosity_info<T> const& info,
        std::vector<T> const& phase_space
    );

    /// Returns the luminosity information.
    luminosity_info<T> const& info() const;

    /// Returns the phase space point.
    std::vector<T> const& phase_space() const;

private:
    luminosity_info<T> info_;
    std::vector<T> phase_space_;
};

/// Creates a phase space generator which generates the points given in the
/// list, in the same ordering. If the list is exhausted, the generator starts
/// with the first point again. Each point has a weight of one.
template <typename T>
std::unique_ptr<phase_space_generator<T>> make_list_phase_space_generator(
    std::vector<list_phase_space_point<T>> const& list_of_phase_space_points
);

}

#endif
