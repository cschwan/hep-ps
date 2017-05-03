#ifndef HEP_PS_MC_OBSERVABLES_HPP
#define HEP_PS_MC_OBSERVABLES_HPP

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
#include "hep/ps/observables.hpp"
#include "hep/ps/mc_distributions.hpp"

#include "hep/mc/multi_channel_point.hpp"
#include "hep/mc/projector.hpp"

#include <memory>
#include <utility>

namespace
{

template <typename T, typename P>
hep::luminosity_info<T> info(hep::multi_channel_point2<T, P> const& x)
{
	return x.map().info();
}

template <typename T, typename P>
hep::luminosity_info<T> info(
	hep::multi_channel_point2<T, std::reference_wrapper<P>> const& x
) {
	return x.map().get().info();
}

}

namespace hep
{

template <typename T>
class mc_observables
{
public:
	mc_observables(
		std::unique_ptr<observables<T>>&& observables,
		initial_state_set set
	)
		: observables_(std::move(observables))
		, distributions_(nullptr)
		, set_(set)
	{
	}

	template <typename P>
	T operator()(
		hep::multi_channel_point2<T, P> const& point,
		hep::projector<T>& projector
	) {
		if (distributions_ == nullptr)
		{
			distributions_ = dynamic_cast <mc_distributions<T>*>
				(&observables_->distributions());
		}

		distributions_->set_projector(projector);
		return observables_->eval(point.coordinates(), info(point), set_);
	}

private:
	std::unique_ptr<observables<T>> observables_;
	mc_distributions<T>* distributions_;
	initial_state_set set_;
};

}

#endif
