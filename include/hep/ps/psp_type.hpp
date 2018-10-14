#ifndef HEP_PS_PSP_TYPE_HPP
#define HEP_PS_PSP_TYPE_HPP

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

/// Enumerator whose two values represent two different phase space points. Each phase space
/// generator generates phase space points in the partonic center-of-mass system, which has to be
/// translated into the corresponding point in the LAB system.
enum class psp_type
{
    pos_rap,
    neg_rap
};

}

#endif
