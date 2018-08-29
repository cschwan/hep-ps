#ifndef HEP_PS_FINAL_STATE_HPP
#define HEP_PS_FINAL_STATE_HPP

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

#include <cstddef>

namespace hep
{

/// Enumeration used to distinguish the final states of matrix elements.
enum class final_state : std::size_t
{
    /// Charged leptons.
    charged_lepton,

    /// Quarks or Gluons.
    quark_gluon,

    /// Photons.
    photon,

    /// Neutrinos.
    neutrino
};

}

#endif
