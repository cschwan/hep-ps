#ifndef HEP_PS_OL_INTEGRATED_MES_HPP
#define HEP_PS_OL_INTEGRATED_MES_HPP

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
#include "hep/ps/dipole_veto.hpp"
#include "hep/ps/final_state.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/photon_dipole_selector.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace hep
{

template <typename T>
class ol_integrated_mes
{
public:
    ol_integrated_mes(
        std::vector<std::string> const& real_processes,
        std::vector<final_state> const& dipole_final_states,
        coupling_order order,
        dipole_veto const& veto = default_dipole_veto(),
        photon_dipole_selector const& selector =
            default_photon_dipole_selector()
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
    std::vector<std::vector<T>> charge_table_;
    std::unordered_multimap<insertion_term, std::tuple<initial_state, int, std::size_t>> mes_;
    std::vector<insertion_term> dipoles_;
    std::vector<final_state> final_states_;
    std::vector<double> ol_m2_;
    std::vector<double> ol_phase_space_;
};

}

#endif
