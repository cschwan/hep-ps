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
#include "hep/ps/particle_type.hpp"

#include <string>
#include <vector>

namespace hep
{

/// Convert the PDG identifier to a \ref final_state.
final_state pdg_id_to_final_state(int id);

/// Converts the PDG identifier to a \ref parton.
parton pdg_id_to_parton(int id);

/// Converts the PDG identifier to a \ref particle_type.
particle_type pdg_id_to_particle_type(int id);

/// Converts the \ref parton `p` to the corresponding particle data group (PDG)
/// ID.
int parton_to_pdg_id(parton p);

/// Checks if the particle with the given PDG identifier carries color.
bool pdg_id_has_color(int id);

/// Returns true if the particle with the given PDF identifier carries an
/// electric charge.
bool pdg_id_has_charge(int id);

/// Returns the PDG identifier of the photon.
int pdg_id_of_photon();

/// Returns the PDG identifier of the gluon.
int pdg_id_of_gluon();

/// Returns true if `id` is the PDG identifier of the gluon.
bool pdg_id_is_gluon(int id);

/// Returns `true` if `id` is the PDG identifier of a quark (or an anti-quark).
bool pdg_id_is_quark(int id);

///
std::pair<initial_state, std::vector<final_state>> pdg_ids_to_states(
	std::vector<int> const& ids
);

/// Converts an OpenLoops process string into a vector of PDG identifiers. This
/// is the inverse function of \ref pdg_ids_to_ol_process_string.
std::vector<int> ol_process_string_to_pdg_ids(std::string const& process);

/// Converts a vector of PDG identifiers to an OpenLoops process string. This
/// is the inverse function of \ref ol_process_string_to_pdg_ids.
std::string pdg_ids_to_ol_process_string(std::vector<int> const& ids);

/// Returns the electromagnetic charge of the particle with the given PDG
/// identifier `id`.
int pdg_id_to_charge_times_three(int id);

}

#endif
