#ifndef HEP_PS_DIAGRAM_HPP
#define HEP_PS_DIAGRAM_HPP

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

#include "hep/ps/propagator.hpp"
#include "hep/ps/vertex.hpp"

#include <vector>

namespace hep
{

/// Class that represents a single Feynman diagram by enumerating its vertices
/// and propagators.
class diagram
{
public:
	/// Constructor.
	diagram(
		std::vector<vertex> const& vertices,
		std::vector<propagator> const& propagators
	);

	/// Returns the vertices of this diagram.
	std::vector<vertex> const& vertices() const;

	/// Returns the propagators of this diagram.
	std::vector<propagator> const& propagators() const;

private:
	std::vector<vertex> vertices_;
	std::vector<propagator> propagators_;
};

}

#endif
