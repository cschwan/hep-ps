#ifndef HEP_PS_ABC_TERMS_HPP
#define HEP_PS_ABC_TERMS_HPP

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

#include "hep/ps/parton.hpp"

namespace hep
{

template <typename T>
struct abc_terms
{
	parton_array<parton_array<T>> a;
	parton_array<parton_array<T>> b;
	parton_array<parton_array<T>> c;
};

}

#endif
