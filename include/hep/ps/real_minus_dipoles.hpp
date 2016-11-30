#ifndef HEP_PS_REAL_MINUS_DIPOLES_HPP
#define HEP_PS_REAL_MINUS_DIPOLES_HPP

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

#include "hep/ps/cut_result.hpp"
#include "hep/ps/initial_state_array.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/particle_type.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

namespace
{

inline std::vector<std::size_t> adjust_indices(
	std::vector<std::size_t> const& indices,
	std::size_t unresolved
) {
	std::vector<std::size_t> result(indices);

	auto const begin = std::find(result.begin(), result.end(), unresolved);
	auto end = result.end();

	assert( begin != end );

	auto next = std::next(begin);

	// decrease all indices following the index for the unresolved by one
	std::transform(next, end, next, [](std::size_t v) { return v - 1; });
	// rotate the unresolved index to the end of the vector
	std::rotate(begin, next, end);
	// remove the unresolved index
	result.erase(std::prev(end));

	return result;
}

inline bool cut_required(hep::initial_state state, hep::cut_result cut)
{
	bool const neg = cut.neg_cutted();
	bool const pos = cut.pos_cutted();

	return (hep::same_initial_states(state) && neg && pos) ||
		(hep::is_negative_ordering(state) && neg) ||
		(hep::is_positive_ordering(state) && pos);
}

}

namespace hep
{

template <typename T, typename M, typename S, typename C, typename R>
class real_minus_dipoles
{
public:
	template <
		typename MatrixElements,
		typename Subtraction,
		typename Cuts,
		typename Recombiner>
	real_minus_dipoles(
		MatrixElements&& matrix_elements,
		Subtraction&& subtraction,
		Cuts&& cuts,
		Recombiner&& recombiner,
		bool inclusive,
		T alpha_min
	)
		: matrix_elements_(std::forward<MatrixElements>(matrix_elements))
		, subtraction_(std::forward<Subtraction>(subtraction))
		, cuts_(std::forward<Cuts>(cuts))
		, recombiner_(std::forward<Recombiner>(recombiner))
		, inclusive_(inclusive)
		, alpha_min_(alpha_min)
	{
	}

	template <typename D1, typename D2>
	initial_state_array<T> operator()(
		std::vector<T> const& real_phase_space,
		T rapidity_shift,
		initial_state_set set,
		D1&& real_differential_info,
		D2&& dipole_differential_info
	) {
		initial_state_array<T> result;
		std::vector<T> aux_phase_phase(real_phase_space.size());

		auto const recombined = recombiner_.recombine(
			real_phase_space,
			aux_phase_phase,
			matrix_elements_.recombination_candidates(),
			1
		);

		bool const is_real = recombined == 1;
		bool const is_inclusive = recombined == 0;

		if (is_real || (inclusive_ && is_inclusive))
		{
			auto const cut_result = cuts_.cut(aux_phase_phase, rapidity_shift,
				is_inclusive);

			if (!cut_result.neg_cutted() || !cut_result.pos_cutted())
			{
				result = matrix_elements_.reals(real_phase_space, set);
				real_differential_info(aux_phase_phase, result);

				if (cut_result.neg_cutted() || cut_result.pos_cutted())
				{
					for (auto const process : set)
					{
						if (cut_required(process, cut_result))
						{
							result.set(process, T());
						}
					}
				}
			}
		}

		aux_phase_phase.resize(real_phase_space.size() - 4);

		// TODO: check for the same dipoles and calculate them only once

		// go through all processes
		for (auto const process : set)
		{
			// go through all dipoles for the current process
			for (auto const dipole : matrix_elements_.dipole_ids(process))
			{
				// map the real phase space on the dipole phase space
				auto const invariants = subtraction_.map_phase_space(
					real_phase_space, aux_phase_phase, dipole);

				if (invariants.adipole < alpha_min_)
				{
					result.set(process, T());
					break;
				}

				auto const dipole_recombination_candidates =
					adjust_indices(matrix_elements_.recombination_candidates(),
						dipole.unresolved());

				auto const dipole_recombined = recombiner_.recombine(
					aux_phase_phase,
					aux_phase_phase,
					dipole_recombination_candidates,
					0
				);

				// check if it passed the recombination
				if (dipole_recombined > 0)
				{
					continue;
				}

				auto const cut_result = cuts_.cut(aux_phase_phase,
					rapidity_shift, false);

				if (cut_required(process, cut_result))
				{
					continue;
				}

				bool const fermion_i = dipole.emitter_type() ==
					particle_type::fermion;
				bool const fermion_j = dipole.unresolved_type() ==
					particle_type::fermion;

				T function;

				if (fermion_i != fermion_j)
				{
					function = subtraction_.fermion_function(dipole,
						invariants);
				}
				else if (fermion_i && fermion_j)
				{
					// TODO: NYI
					assert( false );
				}
				else
				{
					// TODO: NYI
					assert( false );
				}

				T const me = matrix_elements_.dipole(aux_phase_phase, process,
					dipole);
				T const dipole_result = -function * me;

				result.set(process, result.get(process) + dipole_result);
				dipole_differential_info(aux_phase_phase, dipole_result);
			}
		}

		return result;
	}

	initial_state_array<T> operator()(
		std::vector<T> const& real_phase_space,
		T rapidity_shift,
		initial_state_set set
	) {
		return operator()(
			real_phase_space,
			rapidity_shift,
			set,
			[](std::vector<T> const&, initial_state_array<T> const&) {},
			[](std::vector<T> const&, T) {}
		);
	}

	M const& matrix_elements() const
	{
		return matrix_elements_;
	}

	M& matrix_elements()
	{
		return matrix_elements_;
	}

private:
	M matrix_elements_;
	S subtraction_;
	C cuts_;
	R recombiner_;
	bool inclusive_;
	T alpha_min_;
};

template <typename T, typename M, typename S, typename C, typename R>
inline real_minus_dipoles<T, M, S, C, R> make_real_minus_dipoles(
	M&& matrix_elements,
	S&& subtraction,
	C&& cuts,
	R&& recombiner,
	bool inclusive,
	T alpha_min = T()
) {
	return real_minus_dipoles<T, M, S, C, R>(
		std::forward<M>(matrix_elements),
		std::forward<S>(subtraction),
		std::forward<C>(cuts),
		std::forward<R>(recombiner),
		inclusive,
		alpha_min
	);
}

}

#endif
