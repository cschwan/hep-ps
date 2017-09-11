#ifndef HEP_PS_INITIAL_STATE_HPP
#define HEP_PS_INITIAL_STATE_HPP

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

#include "hep/ps/enum.hpp"
#include "hep/ps/parton.hpp"

#include <cassert>

namespace hep
{

/// Enumeration listing all possible partonic initial states. The order matters
/// for the correspoding matrix elements, because of the PDFs. Initial states
/// with transposed partons are automatically taken into account.
HEP_ENUM(initial_state,
	q43_uu, /*!< \f$ Q=4/3 \f$ up/up initial state. */
	q43_cu, /*!< \f$ Q=4/3 \f$ charm/up initial state. */
	q43_cc, /*!< \f$ Q=4/3 \f$ charm/charm initial state. */
	q33_du, /*!< \f$ Q=3/3 \f$ anti-down/up initial state. */
	q33_dc, /*!< \f$ Q=3/3 \f$ anti-down/charm initial state. */
	q33_su, /*!< \f$ Q=3/3 \f$ anti-strange/up initial state. */
	q33_sc, /*!< \f$ Q=3/3 \f$ anti-strange/charm initial state. */
	q23_dd, /*!< \f$ Q=2/3 \f$ anti-down/anti-down initial state. */
	q23_sd, /*!< \f$ Q=2/3 \f$ anti-strange/anti-down initial state. */
	q23_ss, /*!< \f$ Q=2/3 \f$ anti-strange/anti-strange initial state. */
	q23_ug, /*!< \f$ Q=2/3 \f$ up/gluon initial state. */
	q23_cg, /*!< \f$ Q=2/3 \f$ charm/gluon initial state. */
	q13_dg, /*!< \f$ Q=1/3 \f$ anti-down/gluon initial state. */
	q13_sg  /*!< \f$ Q=1/3 \f$ anti-strange/gluon initial state. */
);

HEP_ENUM_ARRAY(initial_state);

HEP_ENUM_SET(initial_state);

// TODO: make function `constexpr` in C++14

/// Returns the second parton of the given initial_state `state`.
inline parton state_parton_one(initial_state state)
{
	switch (state)
	{
	case initial_state::q43_uu:
	case initial_state::q23_ug:
		return parton::up;

	case initial_state::q23_dd:
	case initial_state::q33_du:
	case initial_state::q33_dc:
	case initial_state::q13_dg:
		return parton::anti_down;

	case initial_state::q43_cc:
	case initial_state::q43_cu:
	case initial_state::q23_cg:
		return parton::charm;

	case initial_state::q33_su:
	case initial_state::q23_sd:
	case initial_state::q33_sc:
	case initial_state::q23_ss:
	case initial_state::q13_sg:
		return parton::anti_strange;
	}

	assert( false );
}

// TODO: make function `constexpr` in C++14

/// Returns the second parton of the given initial_state `state`.
inline parton state_parton_two(initial_state state)
{
	switch (state)
	{
	case initial_state::q43_uu:
	case initial_state::q33_du:
	case initial_state::q43_cu:
	case initial_state::q33_su:
		return parton::up;

	case initial_state::q33_dc:
	case initial_state::q43_cc:
	case initial_state::q33_sc:
		return parton::charm;

	case initial_state::q23_dd:
	case initial_state::q23_sd:
		return parton::anti_down;

	case initial_state::q23_ss:
		return parton::anti_strange;

	case initial_state::q23_ug:
	case initial_state::q13_dg:
	case initial_state::q23_cg:
	case initial_state::q13_sg:
		return parton::gluon;
	}

	assert( false );
}

}

#endif
