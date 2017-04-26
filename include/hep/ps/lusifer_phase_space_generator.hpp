#ifndef HEP_PS_LUSIFER_PHASE_SPACE_GENERATOR_HPP
#define HEP_PS_LUSIFER_PHASE_SPACE_GENERATOR_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
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

#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/phase_space_generator.hpp"

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

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_lusifer_phase_space_generator(
	T min_energy,
	T cmf_energy,
	std::vector<std::string> const& processes,
	lusifer_constants<T> const& constants,
	std::size_t extra_random_numbers = 0
);

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_lusifer_phase_space_generator(
	T min_energy,
	T cmf_energy,
	std::string const& process,
	lusifer_constants<T> const& constants,
	std::size_t extra_random_numbers = 0
);

}

#endif
