#ifndef HEP_PS_PHOTON_DIPOLE_SELECTOR_HPP
#define HEP_PS_PHOTON_DIPOLE_SELECTOR_HPP

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

#include <cstddef>
#include <functional>
#include <vector>

namespace hep
{

/// Callback function that allows to select or deselect dipoles used in photon
/// splittings. The first argument is a vector with the PDG ids representing
/// the matrix element that the dipole is proportional to, and the remaining
/// indices are the emitter, unresolved, and spectator index. If the return
/// type is `true`, the corresponding dipole will be included, otherwise it
/// will be discarded.
using photon_dipole_selector = std::function<
    bool(std::vector<int>, std::size_t, std::size_t, std::size_t)>;

/// Selects the first photon dipole which is encountered.
inline photon_dipole_selector default_photon_dipole_selector()
{
    return [](std::vector<int>, std::size_t, std::size_t, std::size_t) {
        return true;
    };
}

}

#endif
