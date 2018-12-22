#ifndef HEP_PS_PS_FUNCTIONS_HPP
#define HEP_PS_PS_FUNCTIONS_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017-2018  Christopher Schwan
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

namespace hep
{

/// Boosts the vector `p` in the direction given by `q`, whose invariant must be `m`. If `inverse`
/// is `true` the inverse boost is applied, otherwise the (uninverted) boost, which t ransforms `p`
/// to the rest frame of `q`.
template <typename T>
void boost(T m, std::array<T, 4>& p, std::array<T, 4> const& q, bool inverse);

/// Rotates the vector `p` first around the 3-axis by `theta` and then around the 2-axis by `phi`.
template <typename T>
void rotate(std::array<T, 4>& p, T phi, T cos_theta);

}

#endif
