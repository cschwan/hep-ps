#ifndef HEP_PS_INDEX_WITH_PARTICLE_CLASS_HPP
#define HEP_PS_INDEX_WITH_PARTICLE_CLASS_HPP

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

#include "hep/ps/particle_class.hpp"

#include <cstddef>

namespace hep
{

/// Class encapsulating a particle index together with a corresponding
/// \ref particle_class.
class index_with_particle_class
{
public:
	/// Constructor.
	index_with_particle_class(
		std::size_t index,
		particle_class particle
	)
		: index_{index}
		, particle_{particle}
	{
	}

	/// Returns the index of the particle.
	std::size_t index() const
	{
		return index_;
	}

	/// Returns the \ref particle_class of the particle with the given
	/// \ref index.
	particle_class particle() const
	{
		return particle_;
	}

private:
	std::size_t index_;
	particle_class particle_;
};

}

#endif
