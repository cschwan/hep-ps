#ifndef HEP_PS_NON_ZERO_DIPOLE_HPP
#define HEP_PS_NON_ZERO_DIPOLE_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017-2018  Christopher Schwan
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

#include "hep/ps/cut_result.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_invariants.hpp"

#include <utility>
#include <vector>

namespace hep
{

/// Class for storing information for the evaluation of a dipole that
/// successfully passed the phase space recombination and cuts.
template <typename T, typename I>
class non_zero_dipole
{
public:
	non_zero_dipole(
		std::size_t index,
		dipole_invariants<T> const& invariants,
		dipole const& dipole,
		cut_result_with_info<I> const& cut_result
	)
		: index_(index)
		, invariants_(invariants)
		, dipole_(dipole)
		, cut_result_(cut_result)
	{
	}

	std::size_t index() const
	{
		return index_;
	}

	dipole_invariants<T> const& invariants() const
	{
		return invariants_;
	}

	hep::dipole const& dipole() const
	{
		return dipole_;
	}

	cut_result_with_info<I> const& cut_result() const
	{
		return cut_result_;
	}

private:
	std::size_t index_;
	dipole_invariants<T> invariants_;
	hep::dipole dipole_;
	cut_result_with_info<I> cut_result_;
};

}

#endif
