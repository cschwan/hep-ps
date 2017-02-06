#ifndef HEP_PS_LUSIFER_PHASE_SPACE_GENERATOR_HPP
#define HEP_PS_LUSIFER_PHASE_SPACE_GENERATOR_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016-2017  Christopher Schwan
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

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace hep
{

template <typename T>
struct lusifer_constants
{
	lusifer_constants(
		T mass_h,
		T width_h,
		T mass_t,
		T width_t,
		T mass_w,
		T width_w,
		T mass_z,
		T width_z
	);

	/// Mass of the Higgs boson.
	T mass_h;

	/// Width of the Higgs boson.
	T width_h;

	/// Mass of the top quark.
	T mass_t;

	/// Width of the top quark.
	T width_t;

	/// Mass of the W bosons.
	T mass_w;

	/// Width of the W boson.
	T width_w;

	/// Mass of the Z boson.
	T mass_z;

	/// Width of the Z boson.
	T width_z;
};

/// A modified version of the phase space generator used in LUSIFER, see
/// \cite Dittmaier:2002ap .
template <typename T>
class lusifer_phase_space_generator
{
public:
	/// The numeric type used to perform phase space generation.
	using numeric_type = T;

	/// Constructor. The parameter `processes` specify the partonic processes
	/// that are used to generate the channels the integrater uses to generate
	/// phase space, with the masses and decay widths specified in `constants`.
	/// The parameter `extra_random_numbers` specifies the number of extra
	/// random numbers that should be generated. This number increases the value
	/// \ref dimensions returns and leaves \ref map_dimensions unchanged. The
	/// random numbers that the generator uses to generate the phase space
	/// points are the first \f$ 3 n - 4 \f$ numbers, where \f$ n \f$ are the
	/// number of final state particles.
	lusifer_phase_space_generator(
		std::vector<std::string> const& processes,
		lusifer_constants<T> const& constants,
		std::size_t extra_random_numbers = 0
	);

	/// Constructor.
	lusifer_phase_space_generator(
		std::string const& process,
		lusifer_constants<T> const& constants,
		std::size_t extra_random_numbers = 0
	);

	/// Move constructor.
	lusifer_phase_space_generator(lusifer_phase_space_generator<T>&& psg);

	/// Destructor.
	~lusifer_phase_space_generator();

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

/// The phase space generator \ref lusifer_phase_space_generator for
/// hadron-hadron collisions.
template <typename T>
using hh_lusifer_phase_space_generator =
	hh_phase_space_generator<lusifer_phase_space_generator<T>>;

}

#endif
