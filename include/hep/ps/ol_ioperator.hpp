#ifndef HEP_PS_OL_IOPERATOR_HPP
#define HEP_PS_OL_IOPERATOR_HPP

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

#include "hep/ps/coupling_order.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/final_state.hpp"
#include "hep/ps/regularization_scheme.hpp"
#include "hep/ps/scales.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace hep
{

///
template <typename T>
class ol_ioperator
{
public:
    ///
    ol_ioperator(
        std::string const& process,
        coupling_order const& order,
        regularization_scheme scheme
    );

    ///
    void alphas(T alphas);

    ///
    std::size_t alphas_power() const;

    ///
    void borns(
        std::vector<T> const& phase_space,
        initial_state_set set,
        std::vector<scales<T>> const& scales,
        std::vector<initial_state_map<T>>& results
    );

    ///
    std::vector<final_state> const& final_states() const;

private:
    int id_;
    coupling_order order_;
    std::vector<final_state> final_states_;
    initial_state state_;
};

}

#endif
