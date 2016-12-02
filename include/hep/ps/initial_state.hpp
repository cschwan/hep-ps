#ifndef HEP_PS_INITIAL_STATE_HPP
#define HEP_PS_INITIAL_STATE_HPP

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

#include <array>
#include <cassert>
#include <cstddef>
#include <initializer_list>

// TODO: make the `enum` an `enum class`

#define HEP_ENUM(name, ...)                                                    \
	enum name : std::size_t                                                    \
	{                                                                          \
		__VA_ARGS__                                                            \
	};                                                                         \
	constexpr std::initializer_list<name> name ## _list()                      \
	{                                                                          \
		return { __VA_ARGS__ };                                                \
	}                                                                          \
	template <typename T>                                                      \
	using name ## _array_ = std::array<T, name ## _list().size()>

namespace hep
{

/// Enumeration listing all possible partonic initial states. The order matters,
/// e.g. there is both \ref q43_uc and \ref q43_cu. If both initial states are
/// the same, e.g. for \ref q43_uu, then there is another value where `x`
/// replaces the `q` in the identifier, e.g. \ref x43_uu. Note that the
/// relative ordering of similar initial states (same partons, but transposed)
/// is important because it determines whether a process belongs to a positive
/// or negative rapidity shift.
HEP_ENUM(initial_state,
	q43_uu,
	x43_uu,
	q43_uc,
	q43_cu,
	q43_cc,
	x43_cc,
	q33_ud,
	q33_du,
	q33_cd,
	q33_dc,
	q33_us,
	q33_su,
	q33_cs,
	q33_sc,
	q23_dd,
	x23_dd,
	q23_ds,
	q23_sd,
	q23_ss,
	x23_ss,
	q23_ug,
	q23_gu,
	q23_cg,
	q23_gc,
	q13_dg,
	q13_gd,
	q13_sg,
	q13_gs
);

// TODO: make function `constexpr` in C++14
inline initial_state swap_initial_state(initial_state state)
{
	switch (state)
	{
	case initial_state::q43_uu: return initial_state::x43_uu;
	case initial_state::q43_cc: return initial_state::x43_cc;
	case initial_state::q23_dd: return initial_state::x23_dd;
	case initial_state::q23_ss: return initial_state::x23_ss;
	case initial_state::x43_uu: return initial_state::q43_uu;
	case initial_state::x43_cc: return initial_state::q43_cc;
	case initial_state::x23_dd: return initial_state::q23_dd;
	case initial_state::x23_ss: return initial_state::q23_ss;

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
	return state < swap_initial_state(state);
}

// TODO: make function `constexpr` in C++14
inline bool state_has_pos_shift(initial_state state)
{
	return state > swap_initial_state(state);
}

}

#endif
