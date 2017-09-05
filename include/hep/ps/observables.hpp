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

#include "hep/mc/multi_channel_point.hpp"
#include "hep/mc/projector.hpp"

#include "hep/ps/initial_state.hpp"
#include "hep/ps/luminosity_info.hpp"

#include <functional>
#include <vector>

namespace
{

template <typename T, typename P>
inline hep::luminosity_info<T> info(hep::multi_channel_point2<T, P> const& x)
{
	return x.map().info();
}

template <typename T, typename P>
inline hep::luminosity_info<T> info(
	hep::multi_channel_point2<T, std::reference_wrapper<P>> const& x
) {
	return x.map().get().info();
}

}

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
	/// `extra_random_numbers`.
	virtual T eval(
		std::vector<T> const& phase_space,
		luminosity_info<T> const& info,
		hep::projector<T>& projector
	) = 0;

	/// Interface for the `hep-mc` Monte Carlo integration routines.
	template <typename P>
	T operator()(
		hep::multi_channel_point2<T, P> const& point,
		hep::projector<T>& projector
	) {
		return eval(point.coordinates(), info(point), projector);
	}
};

}

#endif
