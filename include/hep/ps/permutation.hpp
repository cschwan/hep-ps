#ifndef HEP_PS_PERMUTATION_HPP
#define HEP_PS_PERMUTATION_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017  Christopher Schwan
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
#include <cstddef>

namespace hep
{

template <typename T, std::size_t N>
inline std::array<T, N> inverse_permutation(std::array<T, N> permutation)
{
    // TODO: check if `permutation` is actually a permutation?

    std::array<T, N> inverse;

    for (std::size_t i = 0; i != N; ++i)
    {
        inverse.at(permutation[i]) = i;
    }

    return inverse;
}

}

#endif
