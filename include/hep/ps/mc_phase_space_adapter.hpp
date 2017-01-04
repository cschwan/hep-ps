#ifndef HEP_PS_MC_PHASE_SPACE_ADAPTOR_HPP
#define HEP_PS_MC_PHASE_SPACE_ADAPTOR_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
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

#include "hep/mc/multi_channel_map.hpp"
#include "hep/ps/luminosity_info.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace hep
{

/// Class that bridges the phase space generators from this project with the
/// interface needed by Monte Carlo integrators from `hep-mc`.
template <typename PhaseSpaceGenerator>
class mc_phase_space_adapter
{
private:
	using T = typename PhaseSpaceGenerator::numeric_type;

	PhaseSpaceGenerator psg;
	T cmf_energy_;

public:
	/// Constructor.
	template <typename... Args>
	mc_phase_space_adapter(T cmf_energy, Args&&... args)
		: psg(std::forward<Args>(args)...)
		, cmf_energy_(cmf_energy)
	{
	}

	/// Generates the phase space if `action` is `calculate_momenta` and
	/// calculates the densities otherwise.
	template <typename T>
	T operator()(
		std::size_t channel,
		std::vector<T> const& random_numbers,
		std::vector<T>& momenta,
		std::vector<T>& densities,
		hep::multi_channel_map action
	) {
		if (action == hep::multi_channel_map::calculate_densities)
		{
			return psg.densities(densities);
		}

		psg.generate(random_numbers, momenta, cmf_energy_, channel);

		// value gets ignored
		return T(1.0);
	}

	/// Returns the numbers of channels of `PhaseSpaceGenerator`.
	std::size_t channels() const
	{
		return psg.channels();
	}

	/// Returns the number of dimensions of `PhaseSpaceGenerator`.
	std::size_t dimensions() const
	{
		return psg.dimensions();
	}

	luminosity_info<T> info() const
	{
		return psg.info();
	}

	/// Returns the value of `map_dimensions()` of `PhaseSpaceGenerator`.
	std::size_t map_dimensions() const
	{
		return psg.map_dimensions();
	}
};

}

#endif
