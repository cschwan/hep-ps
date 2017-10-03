#ifndef HEP_PS_PHASE_SPACE_GENERATOR_HPP
#define HEP_PS_PHASE_SPACE_GENERATOR_HPP

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

#include "hep/mc/multi_channel_map.hpp"

#include "hep/ps/luminosity_info.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// Abstract base class of all phase space generators.
template <typename T>
class phase_space_generator
{
public:
	/// The numeric type used to perform phase space generation.
	using numeric_type = T;

	/// Virtual destructor.
	virtual ~phase_space_generator();

	/// Returns the number of channels this generators supports.
	virtual std::size_t channels() const = 0;

	/// Evaluates all densities for the previously generated phase space point
	/// and writes them into `densities`. This vector must have the size given
	/// by \ref channels. The return value is an additional jacobian.
	virtual T densities(std::vector<T>& densities) = 0;

	/// Returns how many random numbers are needed to construct a phase space
	/// point.
	virtual std::size_t dimensions() const = 0;

	/// Uses `random_numbers` to construct a phase space point which is written
	/// as a set of four-vectors into `momenta` using the specified `channel`
	/// for the given center-of-mass frame energy `cmf_energy`.
	virtual void generate(
		std::vector<T> const& random_numbers,
		std::vector<T>& momenta,
		std::size_t channel
	) = 0;

	/// Returns an object of class \ref luminosity_info that contains the energy
	/// of the previously generated phase space point.
	virtual luminosity_info<T> info() const = 0;

	/// Returns the size of the vector `momenta` that has to be passed to
	/// \ref generate.
	virtual std::size_t map_dimensions() const = 0;

	/// Interface for the `hep-mc` Monte Carlo integration routines.
	T operator()(
		std::size_t channel,
		std::vector<T> const& random_numbers,
		std::vector<T>& momenta,
		std::vector<T>& densities,
		hep::multi_channel_map action
	);
};

}

#endif
