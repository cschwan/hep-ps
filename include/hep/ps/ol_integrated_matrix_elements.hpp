#ifndef HEP_PS_OL_INTEGRATED_MATRIX_ELEMENTS_HPP
#define HEP_PS_OL_INTEGRATED_MATRIX_ELEMENTS_HPP

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

#include "hep/ps/correction_type.hpp"
#include "hep/ps/coupling_order.hpp"
#include "hep/ps/final_state.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/insertion_term.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace hep
{

template <typename T>
class ol_integrated_matrix_elements
{
public:
    ol_integrated_matrix_elements(
        std::vector<std::string> const& processes,
        coupling_order const& order,
        correction_type type
    );

    void alphas(T alphas);

    std::size_t alphas_power() const;

    void correlated_me(
        std::vector<T> const& phase_space,
        initial_state_set set,
        std::vector<initial_state_map<T>>& results
    );

    std::vector<final_state> const& final_states() const;

    std::vector<insertion_term> const& insertion_terms() const;

private:
    std::size_t alphas_power_;
    std::unordered_map<int, std::vector<T>> charge_table_;
    std::vector<final_state> final_states_;
#if __GNUC__ > 5
    std::unordered_multimap<initial_state, int> ids_;
#else
    std::unordered_multimap<std::size_t, int> ids_;
#endif
    std::vector<double> ol_m2cc_;
    std::vector<double> ol_phase_space_;
    std::vector<insertion_term> terms_;
    correction_type type_;
};

}

#endif
