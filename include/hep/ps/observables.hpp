#ifndef HEP_PS_OBSERVABLES_HPP
#define HEP_PS_OBSERVABLES_HPP

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

#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/random_numbers.hpp"

#include <vector>

namespace hep
{

/// Abstract base class for all classes that perform the calculation of (parts)
/// of observables.
template <typename T>
class observables
{
public:
	/// Evaluates the observables this instances represents for the given
	/// `phase_space` point and luminosity information in `info`. If the
	/// observables need extra random numbers these will be supplied in
	/// `extra_random_numbers`. If only restricted set of initial states are
	/// needed those will be given in `set`.
	virtual T eval(
		std::vector<T> const& phase_space,
		luminosity_info<T> const& info,
		random_numbers<T>& extra_random_numbers,
		initial_state_set set
	) = 0;
};

}

#endif
