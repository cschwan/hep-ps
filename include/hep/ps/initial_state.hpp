#ifndef HEP_PS_INITIAL_STATE_HPP
#define HEP_PS_INITIAL_STATE_HPP

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

#include "hep/ps/enum.hpp"
#include "hep/ps/parton.hpp"

#include <cassert>
#include <cstddef>

namespace hep
{

/// Enumeration listing all possible partonic initial states. The order matters,
/// e.g. there is both `q43_uc` and `q43_cu`. Note that the relative ordering of
/// similar initial states (same partons, but transposed) is important because
/// it determines whether a process belongs to a positive or negative rapidity
/// shift.
HEP_ENUM(initial_state,
	q43_uu,
	q43_uc,
	q43_cu,
	q43_cc,
	q33_ud,
	q33_du,
	q33_cd,
	q33_dc,
	q33_us,
	q33_su,
	q33_cs,
	q33_sc,
	q23_dd,
	q23_ds,
	q23_sd,
	q23_ss,
	q23_gu,
	q23_ug,
	q23_gc,
	q23_cg,
	q13_gd,
	q13_dg,
	q13_gs,
	q13_sg
);

// TODO: make function `constexpr` in C++14
inline initial_state swap_initial_state(initial_state state)
{
	switch (state)
	{
	case initial_state::q43_uu:
	case initial_state::q43_cc:
	case initial_state::q23_dd:
	case initial_state::q23_ss:
		return state;

	case initial_state::q43_cu: return initial_state::q43_uc;
	case initial_state::q43_uc: return initial_state::q43_cu;
	case initial_state::q33_du: return initial_state::q33_ud;
	case initial_state::q33_ud: return initial_state::q33_du;
	case initial_state::q33_dc: return initial_state::q33_cd;
	case initial_state::q33_cd: return initial_state::q33_dc;
	case initial_state::q33_su: return initial_state::q33_us;
	case initial_state::q33_us: return initial_state::q33_su;
	case initial_state::q33_sc: return initial_state::q33_cs;
	case initial_state::q33_cs: return initial_state::q33_sc;
	case initial_state::q23_sd: return initial_state::q23_ds;
	case initial_state::q23_ds: return initial_state::q23_sd;
	case initial_state::q23_ug: return initial_state::q23_gu;
	case initial_state::q23_gu: return initial_state::q23_ug;
	case initial_state::q23_cg: return initial_state::q23_gc;
	case initial_state::q23_gc: return initial_state::q23_cg;
	case initial_state::q13_dg: return initial_state::q13_gd;
	case initial_state::q13_gd: return initial_state::q13_dg;
	case initial_state::q13_sg: return initial_state::q13_gs;
	case initial_state::q13_gs: return initial_state::q13_sg;
	}

	// implementation error
	assert( false );
}

constexpr bool operator<(initial_state a, initial_state b)
{
	return static_cast <std::size_t> (a) < static_cast <std::size_t> (b);
}

constexpr bool operator>(initial_state a, initial_state b)
{
	return static_cast <std::size_t> (a) > static_cast <std::size_t> (b);
}

// TODO: make function `constexpr` in C++14
inline bool state_has_neg_shift(initial_state state)
{
	return state <= swap_initial_state(state);
}

// TODO: make function `constexpr` in C++14
inline bool state_has_pos_shift(initial_state state)
{
	return state >= swap_initial_state(state);
}

// TODO: make function `constexpr` in C++14
inline parton state_parton_one(initial_state state)
{
	switch (state)
	{
	case initial_state::q43_uu: return parton::up;
	case initial_state::q43_cc: return parton::charm;
	case initial_state::q23_dd: return parton::anti_down;
	case initial_state::q23_ss: return parton::anti_strange;
	case initial_state::q43_cu: return parton::charm;
	case initial_state::q43_uc: return parton::up;
	case initial_state::q33_du: return parton::anti_down;
	case initial_state::q33_ud: return parton::up;
	case initial_state::q33_dc: return parton::anti_down;
	case initial_state::q33_cd: return parton::charm;
	case initial_state::q33_su: return parton::anti_strange;
	case initial_state::q33_us: return parton::up;
	case initial_state::q33_sc: return parton::anti_strange;
	case initial_state::q33_cs: return parton::charm;
	case initial_state::q23_sd: return parton::anti_strange;
	case initial_state::q23_ds: return parton::anti_down;
	case initial_state::q23_ug: return parton::up;
	case initial_state::q23_gu: return parton::gluon;
	case initial_state::q23_cg: return parton::charm;
	case initial_state::q23_gc: return parton::gluon;
	case initial_state::q13_dg: return parton::anti_down;
	case initial_state::q13_gd: return parton::gluon;
	case initial_state::q13_sg: return parton::anti_strange;
	case initial_state::q13_gs: return parton::gluon;
	}

	assert( false );
}

// TODO: make function `constexpr` in C++14
inline parton state_parton_two(initial_state state)
{
	return state_parton_one(swap_initial_state(state));
}

}

#endif
