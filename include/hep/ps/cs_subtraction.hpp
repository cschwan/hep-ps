#ifndef HEP_PS_CS_SUBTRACTION_HPP
#define HEP_PS_CS_SUBTRACTION_HPP

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

#include "hep/ps/ab_terms.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_invariants.hpp"
#include "hep/ps/factorization_scheme.hpp"
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/insertion_term_mode.hpp"
#include "hep/ps/regularization_scheme.hpp"
#include "hep/ps/scales.hpp"
#include "hep/ps/spin_correlation_matrix.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace hep
{

/// Class implementing the Catani-Seymour dipole subtraction for massless
/// particles.
template <typename T>
class cs_subtraction
{
public:
	/// Constructor. Sets the number of colors `nc`, the normalization of the
	/// trace of two generators `tf`, the number of flavors `nf`, the
	/// factorization scheme `fscheme`, and the regularization scheme `rscheme`.
	cs_subtraction(
		T nc,
		T tf,
		std::size_t nf,
		factorization_scheme fscheme,
		regularization_scheme rscheme
	);

	/// Maps the `real_phase_space` onto the `born_phase_space` using the maps
	/// defined for the dipole with `dipole_info`.
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
	T fermion_function(
		dipole const& dipole_info,
		dipole_invariants<T> const& invariants
	);

	/// Returns the finite part of the integrated dipoles together with the
	/// collinear counterterm.
	void insertion_terms(
		insertion_term const& term,
		std::vector<scales<T>> const& scales,
		std::vector<T> const& phase_space,
		T x,
		T eta,
		std::vector<ab_terms<T>>& results
	) const;

	/// Returns the finite part of the integrated dipoles not covered by the
	/// function \ref insertion_terms.
	void insertion_terms2(
		insertion_term const& term,
		std::vector<scales<T>> const& scales,
		std::vector<T> const& phase_space,
		std::vector<T>& results
	) const;

private:
	T nc_;
	T tf_;
	std::size_t nf_;
	factorization_scheme fscheme_;
	regularization_scheme rscheme_;
};

}

#endif
