#ifndef HEP_P_TYPE_JET_ALGORITHM_HPP
#define HEP_P_TYPE_JET_ALGORITHM_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2016-2017  Christopher Schwan
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

#include <cstddef>
#include <vector>

namespace hep
{

/// Implements the p-type jet algorithms (\f$ k_\mathrm{T} \f$, anti-\f$
/// k_\mathrm{T} \f$, and Cambridge-Aachen) described in \cite Cacciari:2008gp .
template <typename T>
class p_type_jet_algorithm
{
public:
	/// Constructor. The choice `p = 1` corresponds to the \f$ k_\mathrm{T}
	/// \f$-algorithm, `p = 0` to the Cambridge-Aachen algorithm, and `p = -1`
	/// to the anti-\f$ k_\mathrm{T} \f$ algorithm.
	p_type_jet_algorithm(T p, T radius);

	/// Use the given `phase_space` and runs the E-scheme recombination
	/// algorithm and writes the result to `recombined_phase_space`. The number
	/// of recombinations is the return value. The phase space can have
	/// arbitrary momenta and the ones that are subject to the recombination
	/// must given as indices in `recombination_candidates`. The value
	/// `max_recombinations` can be used to short-cut the recombination
	/// algorithm. If more than `max_recombinations` are neccessary the return
	/// value is `max_recombinations + 1` and the content of
	/// `recombined_phase_space` is undefined.
	std::size_t recombine(
		std::vector<T> const& phase_space,
		std::vector<T>& recombined_phase_space,
		std::vector<std::size_t> const& recombination_candidates,
		std::size_t max_recombinations
	);

private:
	bool find_jet(std::vector<T>&);

	T p_;
	T radius2_;

	std::vector<std::size_t> candidates_;
	std::vector<T> dib_;
	std::vector<T> dij_;
};

}

#endif
