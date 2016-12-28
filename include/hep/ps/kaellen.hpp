#ifndef HEP_PS_KAELLEN_HPP
#define HEP_PS_KAELLEN_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016  Christopher Schwan
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

#include <type_traits>

namespace hep
{

/// Calculates the Källén function, also known as triangle function, that
/// arises in three-point one-loop functions or particle decays.
template <typename T1, typename T2, typename T3>
inline typename std::common_type<T1, T2, T3>::type kaellen(T1 x, T2 y, T3 z)
{
	return (x - y - z) * (x - y - z) - T1(4.0) * y * z;
}

}

#endif
