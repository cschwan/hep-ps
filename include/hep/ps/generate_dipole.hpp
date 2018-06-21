#ifndef HEP_PS_GENERATE_DIPOLE_HPP
#define HEP_PS_GENERATE_DIPOLE_HPP

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

#include <cstddef>
#include <vector>

namespace hep
{

/// Tries to generate a dipole process for the real process given
/// `process_pdg_ids` at the given coupling order. The parameter `type`
/// determines whether the result should be an EW or a QCD dipole. The integers
/// `i`, `j`, and `k` denote the indices of the emitter-, unresolved, and
/// spectator-particle, respectively. The result is either the PDG identifiers
/// of the dipole process, or an empty vector, in which case the dipole process
/// does not exist.
std::vector<int> generate_dipole(
	std::vector<int> const& process_pdg_ids,
	hep::coupling_order order,
	hep::correction_type type,
	std::size_t i,
	std::size_t j,
	std::size_t k
);

}

#endif
