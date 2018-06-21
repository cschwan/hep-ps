#ifndef HEP_PS_DIPOLE_VETO_HPP
#define HEP_PS_DIPOLE_VETO_HPP

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

#include "hep/ps/final_state.hpp"

#include <functional>
#include <vector>

namespace hep
{

/// Type for dipole vetos. The first argument is a vector with the PDG
/// identifiers of the matrix element that the dipole is proportional to, the
/// second argument is the vector of final states for the dipole.
using dipole_veto = std::function<bool(
	std::vector<int> const&,
	std::vector<final_state> const&
)>;

/// Returns the default dipole veto, which vetos all dipoles that do lead to
/// different final states.
dipole_veto default_dipole_veto();

/// Behaves as \ref default_dipole_veto, but doesn't veto if the process has a
/// \ref final_state::photon instead of a \ref final_state::quark_gluon.
dipole_veto photon_in_jet();

}

#endif
