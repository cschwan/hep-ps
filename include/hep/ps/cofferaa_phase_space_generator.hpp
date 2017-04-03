#ifndef HEP_PS_COFFERAA_PHASE_GENERATOR_HPP
#define HEP_PS_COFFERAA_PHASE_GENERATOR_HPP

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

#include "hep/ps/hh_phase_space_generator.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace hep
{

/// A modified version of the phase space generator used in \f$
/// \mathrm{Coffer}\gamma\gamma \f$, see \cite Bredenstein:2005zk .
template <typename T>
class cofferaa_phase_space_generator
{
public:
	/// The numeric type used to perform phase space generation.
	using numeric_type = T;

	/// Constructor.
	cofferaa_phase_space_generator(
		std::vector<int> const& process,
		lusifer_constants<T> const& constants,
		bool generate_subtraction_channels = true,
		std::size_t extra_random_numbers = 0
	);

	/// Destructor.
	~cofferaa_phase_space_generator();

	/// Returns the number of channels this generators supports.
	std::size_t channels() const;

	/// Evaluates all densities for the previously generated phase space point
	/// and writes them into `densities`. This vector must have the size given
	/// by \ref channels. The return value is an additional jacobian.
	T densities(std::vector<T>& densities);

	/// Returns how many random numbers are needed to construct a phase space
	/// point.
	std::size_t dimensions() const;

	/// Uses `random_numbers` to construct a phase space point which is written
	/// as a set of four-vectors into `momenta` using the specified `channel`
	/// for the given center-of-mass frame energy `cmf_energy`.
	void generate(
		std::vector<T> const& random_numbers,
		std::vector<T>& momenta,
		T cmf_energy,
		std::size_t channel
	);

	/// Returns an object of class \ref luminosity_info that contains the energy
	/// of the previously generated phase space point.
	luminosity_info<T> info() const;

	/// Returns the size of the vector `momenta` that has to be passed to
	/// \ref generate.
	std::size_t map_dimensions() const;

private:
	class impl;
	std::unique_ptr<impl> pimpl;
};

/// The phase space generator \ref cofferaa_phase_space_generator for
/// hadron-hadron collisions.
template <typename T>
using hh_cofferaa_phase_space_generator =
	hh_phase_space_generator<cofferaa_phase_space_generator<T>>;

}

#endif
