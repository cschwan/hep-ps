#ifndef HEP_PS_PHOTON_CONE_RECOMBINER_HPP
#define HEP_PS_PHOTON_CONE_RECOMBINER_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2019  Christopher Schwan
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

#include "hep/ps/final_state.hpp"
#include "hep/ps/recombined_state.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

template <typename T>
class photon_cone_recombiner
{
public:
    photon_cone_recombiner(T radius);

    void recombine(
        std::vector<T> const& phase_space,
        std::vector<final_state> const& final_states,
        std::vector<T>& recombined_phase_space,
        std::vector<recombined_state>& recombined_states
    );

private:
    T radius2_;

    std::vector<std::size_t> candidates_;
    std::vector<T> cones_;
};

}

#endif
