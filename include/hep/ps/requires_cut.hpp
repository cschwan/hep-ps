#ifndef HEP_PS_REQUIRES_CUT_HPP
#define HEP_PS_REQUIRES_CUT_HPP

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

#include "hep/ps/cut_result.hpp"
#include "hep/ps/initial_state.hpp"

namespace hep
{

// TODO: Make `constexpr` in C++14
inline bool requires_cut(hep::initial_state state, hep::cut_result cut)
{
	return (hep::state_has_neg_shift(state) && cut.neg_cutted()) ||
		(hep::state_has_pos_shift(state) && cut.pos_cutted());
}

}

#endif
