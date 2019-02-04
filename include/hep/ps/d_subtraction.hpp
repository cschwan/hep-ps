#ifndef HEP_PS_D_SUBTRACTION_HPP
#define HEP_PS_D_SUBTRACTION_HPP

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

#include "hep/ps/cs_subtraction.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_invariants.hpp"
#include "hep/ps/factorization_scheme.hpp"
#include "hep/ps/regularization_scheme.hpp"
#include "hep/ps/spin_correlation_matrix.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

template <typename T>
class d_subtraction
{
public:
    d_subtraction(
        T nc,
        T tf,
        std::size_t nf,
        factorization_scheme fscheme,
        regularization_scheme rscheme
    );

    bool same_mapping(dipole const& a, dipole const& b) const;

    dipole_invariants<T> map_phase_space(
        std::vector<T> const& real_phase_space,
        std::vector<T>& born_phase_space,
        dipole const& dipole_info
    );

    ///
    spin_correlation_matrix<T> boson_function(
        dipole const& dipole_info,
        dipole_invariants<T> const& invariants,
        std::vector<T> const& phase_space
    );

    ///
    T fermion_function(dipole const& dipole_info, dipole_invariants<T> const& invariants);

private:
    cs_subtraction<T> subtraction_;
};

}

#endif
