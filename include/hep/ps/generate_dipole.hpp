#ifndef HEP_PS_GENERATE_DIPOLE_HPP
#define HEP_PS_GENERATE_DIPOLE_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2018-2019  Christopher Schwan
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
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_split.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// Tries to generate the PDG ids for the dipole process for the real process given by
/// `process_pdg_ids` at the given coupling `order`. The parameter `dipole_info` determines the
/// dipole. The PDG ids are written into `dipole_pdg_ids` and the return value denotes the splitting
/// used for this dipole.
dipole_split generate_dipole(
    std::vector<int> const& process_pdg_ids,
    std::vector<int>& dipole_pdg_ids,
    hep::coupling_order order,
    dipole const& dipole_info
);

}

#endif
