#ifndef HEP_P_TYPE_JET_ALGORITHM_HPP
#define HEP_P_TYPE_JET_ALGORITHM_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2016-2018  Christopher Schwan
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
#include "hep/ps/recombined_state.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// Implements the p-type jet algorithms (\f$ k_\mathrm{T} \f$, anti-\f$ k_\mathrm{T} \f$, and
/// Cambridge-Aachen) described in \cite Cacciari:2008gp .
template <typename T>
class p_type_jet_algorithm
{
public:
    /// Constructor. The choice `p = 1` corresponds to the \f$ k_\mathrm{T} \f$-algorithm, `p = 0`
    /// to the Cambridge-Aachen algorithm, and `p = -1` to the anti-\f$ k_\mathrm{T} \f$ algorithm.
    p_type_jet_algorithm(T p, T radius);

    /// Uses the E-scheme recombination procedure to recombine `phase_space` into
    /// `recombined_phase_space`. Only the momenta of the final states are used that are either a
    /// quark or a gluon.
    void recombine(
        std::vector<T> const& phase_space,
        std::vector<final_state> const& final_states,
        std::vector<T>& recombined_phase_space,
        std::vector<recombined_state>& recombined_states
    );

private:
    T p_;
    T radius2_;

    std::vector<std::size_t> candidates_;
    std::vector<T> dib_;
};

}

#endif
