#ifndef HEP_PS_SILENCE_THIRD_PARTY_LIBS_HPP
#define HEP_PS_SILENCE_THIRD_PARTY_LIBS_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
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

namespace hep
{

/// Prevents third-party libraries, e.g. libraries used for pdf evaluation, from
/// printing their banners. This is useful if multiple processes are started and
/// only one of them should print banners.
void silence_third_party_libs();

}

#endif
