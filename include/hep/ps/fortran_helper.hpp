#ifndef HEP_PS_FORTRAN_HELPER_HPP
#define HEP_PS_FORTRAN_HELPER_HPP

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

#include <vector>

namespace hep
{

/// Assumes `momenta` contains a phase space point in FORTRAN ordering, i.e. first all the energies,
/// then all x-components, and so on. Reorders `momenta` to C++ ordering, i.e. where the momenta are
/// in `0 1 2 3`, `0 1 2 3`, ... ordering.
template <typename T>
void fortran_ordering_to_cpp(std::vector<T>& momenta);

/// Has the opposite effect of \ref fortran_ordering_to_cpp.
template <typename T>
void cpp_ordering_to_fortran(std::vector<T>& momenta);

}

#endif
