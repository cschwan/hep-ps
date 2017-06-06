#ifndef HEP_PS_INV_IDX_HPP
#define HEP_PS_INV_IDX_HPP

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

/// Class representing an invariant in a diagram for a 2 to N scattering
/// process. The indices of the incoming particles must be `1` and `2`, and the
/// remaining ones, i.e. `2`, `3`, etc., are indices of outgoing particles.
class inv_idx
{
public:
	/// Default constructor. Constructs an invariant that is always zero.
	inv_idx();

	/// Constructor.
	inv_idx(std::size_t particle_index, std::size_t n);

	/// Returns a representation of this index.
	std::size_t representation() const;

	/// Returns the number of particles for the process.
	std::size_t n() const;

private:
	std::size_t representation_;
	std::size_t n_;
};

/// Returns `true` if `index` represents an invariant that is s-like, i.e. if
/// it is positive everywhere on phase space. Otherwise this function returns
/// `false` signaling that the invariant is negative everwhere, i.e. a t-like
/// invariant.
bool is_s_like(inv_idx const& index);

/// Returns `true` if `index` represents an invariant that is t-like, i.e. if
/// it is negative everywhere on phase space. Otherwise this function returns
/// `false` signaling that the invariant is positive everwhere, i.e. a s-like
/// invariant.
bool is_t_like(inv_idx const& index);

}

#endif
