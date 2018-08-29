#ifndef HEP_PS_DIPOLE_TYPE_HPP
#define HEP_PS_DIPOLE_TYPE_HPP

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

/// Enumeration for the dipole types.
enum class dipole_type
{
    /// Marks an initial-initial (II) dipole, i.e. a dipole with an initial
    /// state emitter and an initial state spectator.
    initial_initial,

    /// Marks an initial-final (IF) dipole, i.e. a dipole with an initial state
    /// emitter and a final state spectator.
    initial_final,

    /// Marks a final-initial (FI) dipole, i.e. a dipole with a final state
    /// emitter and an initial state spectator.
    final_initial,

    /// Marks a final-final (FF) dipole, i.e. a dipole with a final state
    /// emitter and a final state spectator.
    final_final
};

}

#endif
