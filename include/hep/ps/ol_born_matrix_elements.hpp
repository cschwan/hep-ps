#ifndef HEP_PS_OL_BORN_MATRIX_ELEMENTS_HPP
#define HEP_PS_OL_BORN_MATRIX_ELEMENTS_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2018  Christopher Schwan
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

#include "hep/ps/final_state.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/regularization_scheme.hpp"
#include "hep/ps/scales.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace hep
{

template <typename T>
class ol_born_matrix_elements
{
public:
	ol_born_matrix_elements(
		std::vector<std::string> const& processes,
		std::size_t alphas_power,
		bool loop_mes = false,
		regularization_scheme scheme = regularization_scheme::dim_reg_blha
	);

	void alphas(T alphas);

	std::size_t alphas_power() const;

	void borns(
		std::vector<T> const& phase_space,
		hep::initial_state_set set,
		std::vector<scales<T>> const& scales,
		std::vector<initial_state_map<T>>& results
	);

	std::vector<final_state> const& final_states() const;

private:
	std::size_t alphas_power_;
	bool loop_mes_;
	std::vector<final_state> final_states_;
	std::unordered_multimap<initial_state, int> ids_;
	std::vector<double> ol_phase_space_;
};

}

#endif
