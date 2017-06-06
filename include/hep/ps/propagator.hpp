#ifndef HEP_PS_PROPAGATOR_HPP
#define HEP_PS_PROPAGATOR_HPP

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

#include "hep/ps/inv_idx.hpp"

#include <cstddef>

namespace hep
{

/// Class representing a propagator in a Feynman diagram.
class propagator
{
public:
	/// Constructor.
	propagator(std::size_t type, inv_idx const& invariant);

	/// Returns the type of propagator.
	std::size_t type() const;

	/// Returns an index representation of momenta flowing through this
	/// propagator.
	inv_idx const& invariant() const;

private:
	std::size_t type_;
	inv_idx invariant_;
};

}

#endif
