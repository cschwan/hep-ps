#ifndef HEP_PS_TRIVIAL_CUTTER_HPP
#define HEP_PS_TRIVIAL_CUTTER_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016  Christopher Schwan
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

#include <vector>

namespace hep
{

/// Helper class that never cuts any phase space point. This is useful for
/// testing and calculating total cross sections.
template <typename T>
class trivial_cutter
{
public:
	/// Always returns `false` for positive and negative rapidity shifts.
	cut_result cut(
		std::vector<T> const& phase_space,
		T rapidity_shift,
		bool inclusive_event
	);
};

}

#endif
