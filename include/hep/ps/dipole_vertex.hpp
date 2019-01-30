#ifndef HEP_PS_DIPOLE_VERTEX_HPP
#define HEP_PS_DIPOLE_VERTEX_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2019  Christopher Schwan
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

#include "hep/ps/correction_type.hpp"

namespace hep
{

///
class dipole_vertex
{
public:
    ///
    dipole_vertex(int internal = 0, int external = 0, int unresolved = 0)
        : internal_{internal}
        , external_{external}
        , unresolved_{unresolved}
    {
    }

    ///
    int internal() const
    {
        return internal_;
    }

    ///
    int external() const
    {
        return external_;
    }

    ///
    int unresolved() const
    {
        return unresolved_;
    }

private:
    int internal_;
    int external_;
    int unresolved_;
};

///
bool operator<(dipole_vertex const& a, dipole_vertex const& b);

///
bool operator==(dipole_vertex const& a, dipole_vertex const& b);

///
correction_type correction_type_of(dipole_vertex const& vertex);

}

#endif
