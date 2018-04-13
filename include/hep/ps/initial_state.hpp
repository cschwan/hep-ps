#ifndef HEP_PS_INITIAL_STATE_HPP
#define HEP_PS_INITIAL_STATE_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2016-2018  Christopher Schwan
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

namespace hep
{

/// Enumeration listing all possible partonic initial states. The order matters
/// for the correspoding matrix elements, because of the PDFs. Initial states
/// with transposed partons are automatically taken into account.
HEP_ENUM(initial_state,
	/* Q=4/3 */
	uq_uq, /*!< up up initial state. */
	cq_uq, /*!< charm up initial state. */
	cq_cq, /*!< charm charm initial state. */

	/* Q=1 */
	dx_uq, /*!< anti-down up initial state. */
	dx_cq, /*!< anti-down charm initial state. */
	sx_uq, /*!< anti-strange up initial state. */
	sx_cq, /*!< anti-strange charm initial state. */

	/* Q=2/3 */
	dx_dx, /*!< anti-down anti-down initial state. */
	sx_dx, /*!< anti-strange anti-down initial state. */
	sx_sx, /*!< anti-strange anti-strange initial state. */
	uq_gl, /*!< up gluon initial state. */
	cq_gl, /*!< charm gluon initial state. */

	/* Q=1/3 */
	dx_gl, /*!< anti-down gluon initial state. */
	sx_gl, /*!< anti-strange gluon initial state. */
	uq_dq, /*!< up down initial state. */
	uq_sq, /*!< up strange initial state. */
	cq_dq, /*!< charm down initial state. */
	cq_sq, /*!< charm strange initial state. */

	/* Q=0 */
	uq_ux, /*!< up anti-up initial state. */
	cq_cx, /*!< charm anti-charm initial state. */
	dx_dq, /*!< anti-down down initial state. */
	sx_sq, /*!< anti-strange strange initial state. */
	sx_dq, /*!< anti-strange down initial state. */
	dx_sq, /*!< anti-down strange initial state. */
	cq_ux, /*!< charm anti-up initial state. */
	uq_cx, /*!< up anti-charm initial state. */

	/* Q=-1/3 */
	dx_ux, /*!< anti-down anti-up initial state. */
	sx_ux, /*!< anti-strange anti-up initial state. */
	dx_cx, /*!< anti-down anti-charm initial state. */
	sx_cx  /*!< anti-strange anti-charm initial state. */
);

HEP_ENUM_ARRAY(initial_state);

HEP_ENUM_SET(initial_state);

/// Returns the second parton of the given initial_state `state`.
constexpr parton state_parton_one(initial_state state)
{
	switch (state)
	{
	case uq_uq:
	case uq_gl:
	case uq_dq:
	case uq_sq:
	case uq_ux:
	case uq_cx:
		return parton::up;

	case cq_uq:
	case cq_cq:
	case cq_gl:
	case cq_dq:
	case cq_sq:
	case cq_cx:
	case cq_ux:
		return parton::charm;

	case dx_uq:
	case dx_cq:
	case dx_dx:
	case dx_gl:
	case dx_dq:
	case dx_sq:
	case dx_ux:
	case dx_cx:
		return parton::anti_down;

	case sx_uq:
	case sx_cq:
	case sx_dx:
	case sx_sx:
	case sx_gl:
	case sx_sq:
	case sx_dq:
	case sx_ux:
	case sx_cx:
		return parton::anti_strange;

#if __GNUC__ > 5
	default:
		// if this happens, we didn't cover all the cases
		throw 0;
#endif
	}
}

/// Returns the second parton of the given initial_state `state`.
constexpr parton state_parton_two(initial_state state)
{
	switch (state)
	{
	case uq_uq:
	case cq_uq:
	case dx_uq:
	case sx_uq:
		return parton::up;

	case cq_cq:
	case dx_cq:
	case sx_cq:
		return parton::charm;

	case dx_dx:
	case sx_dx:
		return parton::anti_down;

	case sx_sx:
		return parton::anti_strange;

	case uq_gl:
	case cq_gl:
	case dx_gl:
	case sx_gl:
		return parton::gluon;

	case uq_dq:
	case cq_dq:
	case dx_dq:
	case sx_dq:
		return parton::down;

	case uq_sq:
	case cq_sq:
	case sx_sq:
	case dx_sq:
		return parton::strange;

	case uq_ux:
	case cq_ux:
	case dx_ux:
	case sx_ux:
		return parton::anti_up;

	case cq_cx:
	case uq_cx:
	case dx_cx:
	case sx_cx:
		return parton::anti_charm;

#if __GNUC__ > 5
	default:
		// if this happens, we didn't cover all the cases
		throw 0;
#endif
	}
}

/// Returns a \ref parton_set with all partons that are in `set`.
parton_set partons_in_initial_state_set(initial_state_set set);

}

#endif
