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

#include <array>
#include <bitset>
#include <cstddef>

namespace hep
{

/// Class representing a three- or four-valent vertex in a Feynman diagram.
class vertex
{
public:
	/// Constructor. Creates a vertex connecting particles represented by the
	/// indices `p1`, `p2`, `p3`, and p4`.
	vertex(
		std::size_t p1,
		std::size_t p2,
		std::size_t p3,
		std::size_t p4 = 0
	);

	/// Returns an array of invariants associated with the vertex.
	std::array<std::size_t, 4> const& p() const;

private:
	std::array<std::size_t, 4> p_;
};

}

#endif
