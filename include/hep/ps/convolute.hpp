#ifndef HEP_PS_CONVOLUTE_HPP
#define HEP_PS_CONVOLUTE_HPP

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

#include "hep/ps/cut_result.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/parton.hpp"

namespace hep
{

template <typename T, typename I>
inline neg_pos_results<T> convolute(
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
		auto const sym = (one == two) ? T(0.5) : T(1.0);

		if (!cut.pos_cutted())
		{
			pos += sym * pdfx1[one] * pdfx2[two] * matrix_elements[state];
		}

		if (!cut.neg_cutted())
		{
			neg += sym * pdfx1[two] * pdfx2[one] * matrix_elements[state];
		}
	}

	return { factor * neg , factor * pos };
}

template <typename T, typename I>
inline neg_pos_results<T> convolute(
	parton_array<T> const& pdfx1_neg,
	parton_array<T> const& pdfx2_neg,
	parton_array<T> const& pdfx1_pos,
	parton_array<T> const& pdfx2_pos,
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
		auto const sym = (one == two) ? T(0.5) : T(1.0);

		if (!cut.pos_cutted())
		{
			pos += sym * pdfx1_pos[one] * pdfx2_pos[two] * matrix_elements[state];
		}

		if (!cut.neg_cutted())
		{
			neg += sym * pdfx1_neg[one] * pdfx2_neg[two] * matrix_elements[state];
		}
	}

	return { factor * neg , factor * pos };
}

}

#endif
