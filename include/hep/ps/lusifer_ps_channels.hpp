#ifndef HEP_PS_LUSIFER_PS_CHANNELS_HPP
#define HEP_PS_LUSIFER_PS_CHANNELS_HPP

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

#include "hep/ps/lusifer_constants.hpp"
#include "hep/ps/ps_channel.hpp"

#include <string>
#include <vector>

namespace hep
{

///
template <typename T>
std::vector<ps_channel> lusifer_ps_channels(
    std::vector<std::string> const& processes,
    lusifer_constants<T> const& constants
);

///
template <typename T>
std::vector<ps_channel> lusifer_ps_channels(
    std::string const& process,
    lusifer_constants<T> const& constants
);

}

#endif
