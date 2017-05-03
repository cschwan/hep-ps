#ifndef HEP_PS_MC_PHASE_SPACE_ADAPTOR_HPP
#define HEP_PS_MC_PHASE_SPACE_ADAPTOR_HPP

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

#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/phase_space_generator.hpp"

#include "hep/mc/multi_channel_map.hpp"

#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace hep
{

///
template <typename T>
class mc_phase_space_adapter
{
public:
	/// Constructor.
	mc_phase_space_adapter(std::unique_ptr<phase_space_generator<T>>&& psg)
		: psg_(std::move(psg))
	{
	}

	mc_phase_space_adapter() = default;

	mc_phase_space_adapter(mc_phase_space_adapter<T> const&) = delete;

	mc_phase_space_adapter(mc_phase_space_adapter<T>&& other) = default;

	mc_phase_space_adapter<T>& operator=(mc_phase_space_adapter<T> const&) =
		delete;

	mc_phase_space_adapter<T>& operator=(mc_phase_space_adapter<T>&& other) =
		default;

	/// Generates the phase space if `action` is `calculate_momenta` and
	/// calculates the densities otherwise.
	T operator()(
		std::size_t channel,
		std::vector<T> const& random_numbers,
		std::vector<T>& momenta,
		std::vector<T>& densities,
		hep::multi_channel_map action
	) {
		assert( psg_.get() != nullptr );

		if (action == hep::multi_channel_map::calculate_densities)
		{
			return psg_->densities(densities);
		}

		psg_->generate(random_numbers, momenta, channel);

		// value gets ignored
		return T(1.0);
	}

	std::size_t channels() const
	{
		return psg_->channels();
	}

	std::size_t dimensions() const
	{
		return psg_->dimensions();
	}

	luminosity_info<T> info() const
	{
		return psg_->info();
	}

	std::size_t map_dimensions() const
	{
		return psg_->map_dimensions();
	}

private:
	std::unique_ptr<phase_space_generator<T>> psg_;
};

}

#endif
