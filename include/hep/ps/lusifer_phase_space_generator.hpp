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

#include "hep/ps/lusifer_constants.hpp"
#include "hep/ps/phase_space_generator.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace hep
{

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
