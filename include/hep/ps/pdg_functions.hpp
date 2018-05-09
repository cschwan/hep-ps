#ifndef HEP_PS_PDG_FUNCTIONS_HPP
#define HEP_PS_PDG_FUNCTIONS_HPP

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
#include "hep/ps/parton.hpp"

namespace hep
{

/// Convert the PDG identifier to a \ref final_state.
final_state pdg_id_to_final_state(int id);

/// Converts the PDG identifier to a \ref parton.
parton pdg_id_to_parton(int id);

/// Converts the \ref parton `p` to the corresponding particle data group (PDG)
/// ID.
int parton_to_pdg_id(parton p);

/// Convert two PDG identifiers to an \ref initial_state.
inline initial_state pdg_ids_to_initial_state(int id1, int id2)
{
	return partons_to_initial_state(pdg_id_to_parton(id1),
		pdg_id_to_parton(id2));
}

}

#endif
