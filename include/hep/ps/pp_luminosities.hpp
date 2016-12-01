#ifndef HEP_PS_PP_LUMINOSITIES_HPP
#define HEP_PS_PP_LUMINOSITIES_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016  Christopher Schwan
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

#include "hep/ps/initial_state_array.hpp"

#include <cstddef>
#include <memory>
#include <string>

namespace hep
{

template <typename T>
class pp_luminosities
{
public:
	// TODO: replace `string` with `string_view` in C++17

	pp_luminosities(std::string const& pdf_name, std::size_t pdf_member);

	pp_luminosities(pp_luminosities&& luminosities);

	~pp_luminosities();

	T alphas(T scale);

	initial_state_array<T> pdfs(T x1, T x2, T scale);

private:
	class impl;
	std::unique_ptr<impl> pimpl;
};

}

#endif
