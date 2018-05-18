#ifndef HEP_PS_DIPOLE_WITH_SET_HPP
#define HEP_PS_DIPOLE_WITH_SET_HPP

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

#include "hep/ps/correction_type.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_type.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/particle_type.hpp"

#include <cstddef>

namespace hep
{

/// Class that wraps a \ref dipole together with an \ref initial_state_set.
class dipole_with_set
{
public:
	/// Constructor.
	dipole_with_set(
		std::size_t emitter,
		std::size_t unresolved,
		std::size_t spectator,
		particle_type emitter_type,
		particle_type unresolved_type,
		particle_type spectator_type,
		correction_type corr_type,
		initial_state_set set
	)
		: dipole_(emitter, unresolved, spectator, emitter_type, unresolved_type,
			spectator_type, corr_type)
		, set_(set)
	{
	}

	/// Returns the dipole.
	hep::dipole const& dipole() const
	{
		return dipole_;
	}

	/// Returns the set.
	initial_state_set set() const
	{
		return set_;
	}

private:
	hep::dipole dipole_;
	initial_state_set set_;
};

}

#endif
