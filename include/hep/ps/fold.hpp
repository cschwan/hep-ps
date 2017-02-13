#ifndef HEP_PS_FOLD_HPP
#define HEP_PS_FOLD_HPP

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

#include "hep/ps/initial_state_array.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/neg_pos_results.hpp"

namespace hep
{

template <typename T>
inline T fold(
	initial_state_array<T> const& a,
	initial_state_array<T> const& b,
	initial_state_set set
) {
	T result{};

	for (auto const state : set)
	{
		result += a.get(state) * b.get(state);
	}

	return result;
}

template <typename T>
inline neg_pos_results<T> fold(
	initial_state_array<T> const& a,
	initial_state_array<T> const& b,
	initial_state_set set,
	T factor,
	cut_result cut
) {
	T neg{};
	T pos{};

	for (auto const state : set)
	{
		if (state_has_neg_shift(state) && !cut.neg_cutted())
		{
			neg += a.get(state) * b.get(state);
		}

		if (state_has_pos_shift(state) && !cut.pos_cutted())
		{
			pos += a.get(state) * b.get(state);
		}
	}

	return { factor * neg, factor * pos };
}

template <typename T>
inline neg_pos_results<T> fold(
	initial_state_array<T> const& a,
	T b,
	initial_state state,
	T factor
) {
	T neg{};
	T pos{};

	if (state_has_neg_shift(state))
	{
		neg += a.get(state) * b;
	}
	else
	{
		pos += a.get(state) * b;
	}

	return { factor * neg, factor * pos };
}

}

#endif
