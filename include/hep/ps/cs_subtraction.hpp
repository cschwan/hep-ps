#ifndef HEP_PS_CS_SUBTRACTION_HPP
#define HEP_PS_CS_SUBTRACTION_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
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
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/renormalization_scheme.hpp"
#include "hep/ps/scales.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

template <typename T>
class cs_subtraction
{
public:
	cs_subtraction(
		T n,
		T tf,
		T nf,
		factorization_scheme fscheme,
		renormalization_scheme rscheme
	);

	dipole_invariants<T> map_phase_space(
		std::vector<T> const& real_phase_space,
		std::vector<T>& born_phase_space,
		dipole const& dipole_info
	);

	T fermion_function(
		dipole const& dipole_info,
		dipole_invariants<T> const& invariants
	);

	abc_terms<T> insertion_terms(
		insertion_term const& term,
		scales<T> const& mu,
		std::vector<T> const& phase_space,
		T x,
		T eta,
		std::size_t initial_state
	) const;

	T insertion_terms2(
		insertion_term const& term,
		scales<T> const& mu,
		std::vector<T> const& phase_space
	) const;

private:
	T n_;
	T tf_;
	T nf_;
	factorization_scheme fscheme_;
	renormalization_scheme rscheme_;
};

}

#endif
