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

template <typename T>
class p_type_jet_algorithm
{
public:
	p_type_jet_algorithm(T p, T radius);

	std::size_t recombine(
		std::vector<T> const& phase_space,
		std::vector<T>& recombined_phase_space,
		std::vector<std::size_t> const& recombination_candidates,
		std::size_t max_recombinations
	);

private:
	bool find_jet(std::vector<T>&);

	T p_;
	T radius_;

	std::vector<std::size_t> candidates_;
	std::vector<T> dib_;
	std::vector<T> dij_;
};

}

#endif
