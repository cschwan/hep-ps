#ifndef HEP_PS_LUSIFER_PS_FUNCTIONS_HPP
#define HEP_PS_LUSIFER_PS_FUNCTIONS_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2018  Christopher Schwan
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

namespace hep
{

/// Returns the jacobian for \ref lusifer_invariant_map.
template <typename T>
T lusifer_invariant_jacobian(T power, T mass, T width, T x, T xmin, T xmax);

/// Maps a random number `x` to an interval given by `xmin` and `xmax` according a map which is
/// defined by the three parameters `power`, `mass`, and `width`. The jacobian is computed by the
/// function \ref lusifer_invariant_map.
template <typename T>
T lusifer_invariant_map(T power, T mass, T width, T x, T xmin, T xmax);

}

#endif
