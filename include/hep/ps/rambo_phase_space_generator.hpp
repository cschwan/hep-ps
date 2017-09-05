#ifndef HEP_PS_RAMBO_PHASE_SPACE_GENERATOR_HPP
#define HEP_PS_RAMBO_PHASE_SPACE_GENERATOR_HPP

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

#include "hep/ps/phase_space_generator.hpp"

#include <cstddef>
#include <memory>

namespace hep
{

/// Creates a new generator that generates phase space using the RAMBO
/// algorithm. See \cite Kleiss:1985gy .
template <typename T>
std::unique_ptr<phase_space_generator<T>> make_rambo_phase_space_generator(
	T min_energy,
	T cmf_energy,
	std::size_t final_state_particles,
	std::size_t extra_random_numbers = 0
);

}

#endif
