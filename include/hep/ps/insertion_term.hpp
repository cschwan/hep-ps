#ifndef HEP_PS_INSERTION_TERM_HPP
#define HEP_PS_INSERTION_TERM_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
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

#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/particle_type.hpp"

#include <cstddef>

namespace hep
{

/// Classifies a finite insertion term according to the emitter and spectator
/// particle, similarly to \ref dipole.
class insertion_term
{
public:
	/// Constructor.
	insertion_term(
		std::size_t emitter,
		std::size_t spectator,
		particle_type emitter_type,
		insertion_term_type type
	)
		: emitter_(emitter)
		, spectator_(spectator)
		, emitter_type_(emitter_type)
		, type_(type)
	{
	}

	/// Returns the index of the emitter particle.
	std::size_t emitter() const
	{
		return emitter_;
	}

	/// Returns the index of the spectator particle.
	std::size_t spectator() const
	{
		return spectator_;
	}

	/// Returns the type of the emitter particle.
	particle_type emitter_type() const
	{
		return emitter_type_;
	}

	/// Returns the type of this insertion term.
	insertion_term_type type() const
	{
		return type_;
	}

private:
	std::size_t emitter_;
	std::size_t spectator_;
	particle_type emitter_type_;
	insertion_term_type type_;
};

}

#endif
