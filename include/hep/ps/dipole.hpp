#ifndef HEP_PS_DIPOLE_HPP
#define HEP_PS_DIPOLE_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2016-2017  Christopher Schwan
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

#include "hep/ps/dipole_type.hpp"
#include "hep/ps/particle_type.hpp"

#include <cstddef>

namespace hep
{

class dipole
{
public:
	dipole(
		std::size_t emitter,
		std::size_t unresolved,
		std::size_t spectator,
		particle_type emitter_type,
		particle_type unresolved_type,
		particle_type spectator_type,
		dipole_type type
	)
		: emitter_(emitter)
		, unresolved_(unresolved)
		, spectator_(spectator)
		, emitter_type_(emitter_type)
		, unresolved_type_(unresolved_type)
		, spectator_type_(spectator_type)
		, type_(type)
	{
	}

	std::size_t emitter() const
	{
		return emitter_;
	}

	std::size_t unresolved() const
	{
		return unresolved_;
	}

	std::size_t spectator() const
	{
		return spectator_;
	}

	particle_type emitter_type() const
	{
		return emitter_type_;
	}

	particle_type unresolved_type() const
	{
		return unresolved_type_;
	}

	particle_type spectator_type() const
	{
		return spectator_type_;
	}

	dipole_type type() const
	{
		return type_;
	}

private:
	std::size_t emitter_;
	std::size_t unresolved_;
	std::size_t spectator_;
	particle_type emitter_type_;
	particle_type unresolved_type_;
	particle_type spectator_type_;
	dipole_type type_;
};

bool operator<(dipole const& a, dipole const& b);

bool operator==(dipole const& a, dipole const& b);

}

#endif
