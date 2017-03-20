#ifndef HEP_PS_DIPOLE_INVARIANTS_HPP
#define HEP_PS_DIPOLE_INVARIANTS_HPP

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
struct dipole_invariants
{
	dipole_invariants() = default;

	dipole_invariants(T one, T two, T sij, T adipole)
		: one(one)
		, two(two)
		, sij(sij)
		, adipole(adipole)
	{
	}

	T one;
	T two;
	T sij;
	T adipole;
};

}

#endif
