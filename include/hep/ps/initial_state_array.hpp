#ifndef HEP_PS_INITIAL_STATE_ARRAY_HPP
#define HEP_PS_INITIAL_STATE_ARRAY_HPP

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

#include "hep/ps/initial_state.hpp"

#include <array>
#include <cstddef>

namespace hep
{

template <typename T>
class initial_state_array
{
public:
	initial_state_array()
		// ZERO-initialize the array (!?)
		: array_{{}}
	{
		// TODO: is this really needed?
		array_.fill(T());
	}

	T get(initial_state state) const
	{
		return array_[static_cast <std::size_t> (state)];
	}

	void set(initial_state state, T const& value)
	{
		array_[static_cast <std::size_t> (state)] = value;
	}

private:
	initial_state_array_<T> array_;
};

}

#endif
