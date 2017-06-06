#ifndef HEP_PS_VERTEX_HPP
#define HEP_PS_VERTEX_HPP

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

#include <array>

namespace hep
{

/// Class representing a three- or four-valent vertex in a Feynman diagram.
class vertex
{
public:
	/// Constructor. Creates a vertex connecting particles with invariants
	/// `inv1`, `inv2`, `inv3`, and `inv4`. If one of the invariants represents
	/// the null-invariant this vertex is a three-valent vertex, otherwise a
	/// four-valent vertex.
	vertex(
		inv_idx const& inv1,
		inv_idx const& inv2,
		inv_idx const& inv3,
		inv_idx const& inv4 = inv_idx()
	);

	/// Returns an array of invariants associated with the vertex.
	std::array<inv_idx, 4> const& invariants() const;

private:
	std::array<inv_idx, 4> invariants_;
};

}

#endif
