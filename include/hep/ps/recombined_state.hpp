#ifndef HEP_PS_RECOMBINED_STATE_HPP
#define HEP_PS_RECOMBINED_STATE_HPP

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

/// Enumeration listing the possible objects in a detector.
enum class recombined_state : std::size_t
{
	/// Charged leptons that are possibly recombined with photons.
	dressed_lepton,

	/// Quarks or Gluons or recombinations thereof that form a jet.
	jet,

	/// An isolated photon.
	isolated_photon,

	/// The sum of all neutrino momenta.
	missing_momentum
};

}

#endif
