#ifndef HEP_PS_CONSTANTS_HPP
#define HEP_PS_CONSTANTS_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2016-2017  Christopher Schwan
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

template <typename T>
struct constants
{
	/// The value of \f$ ( \hbar c )^2 \f$ in \f$ \text{GeV}^2 \text{pbarn} \f$
	/// from the 2016 edition of the PDG booklet.
	static constexpr T hbarc2 = T(3.893793656e8);
};

}

#endif
