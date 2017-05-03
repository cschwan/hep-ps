#ifndef HEP_PS_MC_DISTRIBUTIONS_HPP
#define HEP_PS_MC_DISTRIBUTIONS_HPP

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

#include "hep/ps/cut_result.hpp"
#include "hep/ps/distributions.hpp"
#include "hep/ps/event_type.hpp"
#include "hep/ps/neg_pos_results.hpp"

#include "hep/mc/projector.hpp"

namespace hep
{

template <typename T>
class mc_distributions : public distributions<T>
{
public:
	/// Set the projector for the next call
	virtual void set_projector(hep::projector<T>& projector) = 0;
};

template <typename T>
class mc_trivial_distributions : public mc_distributions<T>
{
public:
	void set_projector(hep::projector<T>&) override
	{
	}

	template <typename I>
	void operator()(
		std::vector<T> const&,
		cut_result_with_info<I> const&,
		neg_pos_results<T> const&,
		T,
		event_type
	) {
	}
};

}

#endif