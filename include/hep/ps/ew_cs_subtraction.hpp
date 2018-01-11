#ifndef HEP_PS_EW_CS_SUBTRACTION_HPP
#define HEP_PS_EW_CS_SUBTRACTION_HPP

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

#include "hep/ps/cs_subtraction.hpp"

namespace hep
{

template <typename T>
class ew_cs_subtraction : public cs_subtraction<T>
{
public:
	ew_cs_subtraction(
		T nf,
		factorization_scheme fscheme,
		renormalization_scheme rscheme
	);
};

}

#endif
