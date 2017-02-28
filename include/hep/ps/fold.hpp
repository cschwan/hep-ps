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

#include "hep/ps/cut_result.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/parton.hpp"

namespace hep
{

template <typename T, typename I>
inline neg_pos_results<T> fold(
	parton_array<T> const& pdfx1,
	parton_array<T> const& pdfx2,
	initial_state_array<T> const& matrix_elements,
	initial_state_set set,
	T factor,
	cut_result_with_info<I> const& cut
) {
	T neg{};
	T pos{};

	for (auto const state : set)
	{
		auto const one = state_parton_one(state);
		auto const two = state_parton_two(state);

		if (state_has_neg_shift(state) && !cut.neg_cutted())
		{
			neg += pdfx1[one] * pdfx2[two] * matrix_elements[state];
		}

		if (state_has_pos_shift(state) && !cut.pos_cutted())
		{
			pos += pdfx1[one] * pdfx2[two] * matrix_elements[state];
		}
	}

	return { factor * neg , factor * pos };
}

template <typename T>
inline neg_pos_results<T> fold(
	parton_array<T> const& pdfx1,
	parton_array<T> const& pdfx2,
	T matrix_element,
	initial_state state,
	T factor
) {
	T neg{};
	T pos{};

	auto const one = state_parton_one(state);
	auto const two = state_parton_two(state);

	if (state_has_neg_shift(state))
	{
		neg += pdfx1[one] * pdfx2[two] * matrix_element;
	}

	if (state_has_pos_shift(state))
	{
		pos += pdfx1[one] * pdfx2[two] * matrix_element;
	}

	return { factor * neg , factor * pos };
}

}

#endif