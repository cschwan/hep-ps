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

/// Structure capturing the invariants for a dipole for a specific phase space
/// point.
template <typename T>
struct dipole_invariants
{
    /// Constructor. Initializes this object
    dipole_invariants(T one = T(), T two = T(), T sij = T(), T alpha = T())
        : one{one}
        , two{two}
        , sij{sij}
        , alpha{alpha}
    {
    }

    /// Invariant.
    T one;

    /// Other invariant.
    T two;

    /// Propagator invariant of the dipole.
    T sij;

    /// The value of the parameter \f$ \alpha \f$ on which a technical cut is
    /// applied if the value is too small.
    T alpha;
};

}

#endif
