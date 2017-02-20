#ifndef HEP_PS_CS_SUBTRACTION_HPP
#define HEP_PS_CS_SUBTRACTION_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016-2017  Christopher Schwan
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

#include "hep/ps/abc_terms.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_invariants.hpp"
#include "hep/ps/factorization_scheme.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

template <typename T>
class cs_subtraction
{
public:
	cs_subtraction(T n, T tf, T nf, factorization_scheme scheme);

	dipole_invariants<T> map_phase_space(
		std::vector<T> const& real_phase_space,
		std::vector<T>& born_phase_space,
		dipole const& dipole_info
	);

	T fermion_function(
		dipole const& dipole_info,
		dipole_invariants<T> const& invariants
	);

	///
	abc_terms<T> fermion_abc_born(T x, T eta) const;

	///
	abc_terms<T> fermion_abc_final_initial(
		T x,
		T eta,
		particle_type type_emitter
	) const;

	///
	abc_terms<T> fermion_abc_initial_final(
		T x,
		T eta,
		T mu2,
		std::vector<T> const& phase_space,
		std::size_t emitter,
		std::size_t spectator
	) const;

private:
	T n_;
	T tf_;
	T nf_;
	factorization_scheme scheme_;
};

}

#endif
