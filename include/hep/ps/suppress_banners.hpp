#ifndef HEP_PS_SUPRESS_BANNERS_HPP
#define HEP_PS_SUPRESS_BANNERS_HPP

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

namespace hep
{

/// Prevents all third-party libraries used by `hep-ps` from printing their
/// banners. This is useful if many processes are started and only one of them
/// should print banners.
void suppress_banners(bool suppress);

/// Returns `true` if third-party libraries should be silenced.
bool& suppress_banners();

}

#endif
