#ifndef HEP_PS_PS_HPP
#define HEP_PS_PS_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016  Christopher Schwan
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

extern "C"
{

/// C-function that lets build systems find the `hep-ps` library. This function
/// does nothing.
void find_hep_ps();

}

#include "hep/ps/constants.hpp"
#include "hep/ps/cs_subtraction.hpp"
#include "hep/ps/cut_result.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_invariants.hpp"
#include "hep/ps/dipole_type.hpp"
#include "hep/ps/fortran_helper.hpp"
#include "hep/ps/hh_phase_space_generator.hpp"
#include "hep/ps/initial_state_array.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/kaellen.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "hep/ps/no_cutter.hpp"
#include "hep/ps/no_recombiner.hpp"
#include "hep/ps/observables_born_like.hpp"
#include "hep/ps/observables_real.hpp"
#include "hep/ps/particle_type.hpp"
#include "hep/ps/pp_luminosities.hpp"
#include "hep/ps/p_type_jet_algorithm.hpp"
#include "hep/ps/requires_cut.hpp"
#include "hep/ps/scales.hpp"

#endif
