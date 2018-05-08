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
	cq_cq, /*!< \f$ Q=4/3 \f$ charm charm initial state. */
	cq_cx, /*!< \f$ Q=0 \f$ charm anti-charm initial state. */
	cq_dq, /*!< \f$ Q=1/3 \f$ charm down initial state. */
	cq_dx, /*!< \f$ Q=1 \f$ charm anti-down initial state. */
	cq_gl, /*!< \f$ Q=2/3 \f$ charm gluon initial state. */
	cq_ph, /*!< \f$ Q=2/3 \f$ charm photon initial state. */
	cq_sq, /*!< \f$ Q=1/3 \f$ charm strange initial state. */
	cq_sx, /*!< \f$ Q=1 \f$ charm anti-strange initial state. */
	cq_uq, /*!< \f$ Q=4/3 \f$ charm up initial state. */
	cq_ux, /*!< \f$ Q=0 \f$ charm anti-up initial state. */
	cx_cq, /*!< \f$ Q=0 \f$ anti-charm charm initial state. */
	cx_dx, /*!< \f$ Q=-1/3 \f$ anti-charm anti-down initial state. */
	cx_sx, /*!< \f$ Q=-1/3 \f$ anti-charm anti-strange initial state. */
	cx_uq, /*!< \f$ Q=0 \f$ anti-charm up initial state. */
	dq_cq, /*!< \f$ Q=1/3 \f$ down charm initial state. */
	dq_dx, /*!< \f$ Q=0 \f$ down anti-down initial state. */
	dq_sx, /*!< \f$ Q=0 \f$ down anti-strange initial state. */
	dq_uq, /*!< \f$ Q=1/3 \f$ down up initial state. */
	dx_cq, /*!< \f$ Q=1 \f$ anti-down charm initial state. */
	dx_cx, /*!< \f$ Q=-1/3 \f$ anti-down anti-charm initial state. */
	dx_dq, /*!< \f$ Q=0 \f$ anti-down down initial state. */
	dx_dx, /*!< \f$ Q=2/3 \f$ anti-down anti-down initial state. */
	dx_gl, /*!< \f$ Q=1/3 \f$ anti-down gluon initial state. */
	dx_sq, /*!< \f$ Q=0 \f$ anti-down strange initial state. */
	dx_sx, /*!< \f$ Q=2/3 \f$ anti-strange anti-down initial state. */
	dx_uq, /*!< \f$ Q=1 \f$ anti-down up initial state. */
	dx_ux, /*!< \f$ Q=-1/3 \f$ anti-down anti-up initial state. */
	gl_cq, /*!< \f$ Q=2/3 \f$ gluon charm initial state. */
	gl_dx, /*!< \f$ Q=1/3 \f$ gluon anti-down initial state. */
	gl_gl, /*!< \f$ Q=0 \f$ gluon gluon initial state. */
	gl_sx, /*!< \f$ Q=1/3 \f$ gluon anti-strange initial state. */
	gl_uq, /*!< \f$ Q=2/3 \f$ gluon up initial state. */
	ph_cq, /*!< \f$ Q=2/3 \f$ photon charm initial state. */
	ph_ph, /*!< \f$ Q=0 \f$ photon photon initial state. */
	ph_uq, /*!< \f$ Q=2/3 \f$ photon up initial state. */
	sq_cq, /*!< \f$ Q=1/3 \f$ strange charm initial state. */
	sq_dx, /*!< \f$ Q=0 \f$ strange anti-down initial state. */
	sq_sx, /*!< \f$ Q=0 \f$ strange anti-strange initial state. */
	sq_uq, /*!< \f$ Q=1/3 \f$ strange up initial state. */
	sx_cq, /*!< \f$ Q=1 \f$ anti-strange charm initial state. */
	sx_cx, /*!< \f$ Q=-1/3 \f$ anti-strange anti-charm initial state. */
	sx_dq, /*!< \f$ Q=0 \f$ anti-strange down initial state. */
	sx_dx, /*!< \f$ Q=2/3 \f$ anti-strange anti-down initial state. */
	sx_gl, /*!< \f$ Q=1/3 \f$ anti-strange gluon initial state. */
	sx_sq, /*!< \f$ Q=0 \f$ anti-strange strange initial state. */
	sx_sx, /*!< \f$ Q=2/3 \f$ anti-strange anti-strange initial state. */
	sx_uq, /*!< \f$ Q=1 \f$ anti-strange up initial state. */
	sx_ux, /*!< \f$ Q=-1/3 \f$ anti-strange anti-up initial state. */
	uq_cq, /*!< \f$ Q=4/3 \f$ up charm initial state. */
	uq_cx, /*!< \f$ Q=0 \f$ up anti-charm initial state. */
	uq_dq, /*!< \f$ Q=1/3 \f$ up down initial state. */
	uq_dx, /*!< \f$ Q=1 \f$ up anti-down initial state. */
	uq_gl, /*!< \f$ Q=2/3 \f$ up gluon initial state. */
	uq_ph, /*!< \f$ Q=2/3 \f$ up photon initial state. */
	uq_sq, /*!< \f$ Q=1/3 \f$ up strange initial state. */
	uq_sx, /*!< \f$ Q=1 \f$ up anti-strange initial state. */
	uq_uq, /*!< \f$ Q=4/3 \f$ up up initial state. */
	uq_ux, /*!< \f$ Q=0 \f$ up anti-up initial state. */
	ux_cq, /*!< \f$ Q=0 \f$ anti-up charm initial state. */
	ux_dx, /*!< \f$ Q=-1/3 \f$ anti-up anti-down initial state. */
	ux_sx, /*!< \f$ Q=-1/3 \f$ anti-up anti-strange initial state. */
	ux_uq  /*!< \f$ Q=0 \f$ anti-up up initial state. */
);

HEP_ENUM_MAP(initial_state);

HEP_ENUM_SET(initial_state);

/// Returns the second parton of the given initial_state `state`.
constexpr parton state_parton_one(initial_state state)
{
	switch (state)
	{
	case cq_cq:
	case cq_cx:
	case cq_dq:
	case cq_dx:
	case cq_gl:
	case cq_ph:
	case cq_sq:
	case cq_sx:
	case cq_uq:
	case cq_ux:
		return parton::charm;

	case cx_cq:
	case cx_dx:
	case cx_sx:
	case cx_uq:
		return parton::anti_charm;

	case dq_cq:
	case dq_dx:
	case dq_sx:
	case dq_uq:
		return parton::down;

	case dx_cq:
	case dx_cx:
	case dx_dq:
	case dx_dx:
	case dx_gl:
	case dx_sq:
	case dx_sx:
	case dx_uq:
	case dx_ux:
		return parton::anti_down;

	case gl_cq:
	case gl_dx:
	case gl_gl:
	case gl_sx:
	case gl_uq:
		return parton::gluon;

	case ph_cq:
	case ph_ph:
	case ph_uq:
		return parton::photon;

	case sq_cq:
	case sq_dx:
	case sq_sx:
	case sq_uq:
		return parton::strange;

	case sx_cq:
	case sx_cx:
	case sx_dq:
	case sx_dx:
	case sx_gl:
	case sx_sq:
	case sx_sx:
	case sx_uq:
	case sx_ux:
		return parton::anti_strange;

	case uq_cq:
	case uq_cx:
	case uq_dq:
	case uq_dx:
	case uq_gl:
	case uq_ph:
	case uq_sq:
	case uq_sx:
	case uq_uq:
	case uq_ux:
		return parton::up;

	case ux_cq:
	case ux_dx:
	case ux_sx:
	case ux_uq:
		return parton::anti_up;

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
	case cq_cq:
	case cx_cq:
	case dq_cq:
	case dx_cq:
	case gl_cq:
	case ph_cq:
	case sq_cq:
	case sx_cq:
	case uq_cq:
	case ux_cq:
		return parton::charm;

	case cq_cx:
	case dx_cx:
	case sx_cx:
	case uq_cx:
		return parton::anti_charm;

	case cq_dq:
	case dx_dq:
	case sx_dq:
	case uq_dq:
		return parton::down;

	case cq_dx:
	case cx_dx:
	case dq_dx:
	case dx_dx:
	case gl_dx:
	case sq_dx:
	case sx_dx:
	case uq_dx:
	case ux_dx:
		return parton::anti_down;

	case cq_gl:
	case dx_gl:
	case gl_gl:
	case sx_gl:
	case uq_gl:
		return parton::gluon;

	case cq_ph:
	case ph_ph:
	case uq_ph:
		return parton::photon;

	case cq_sq:
	case dx_sq:
	case sx_sq:
	case uq_sq:
		return parton::strange;

	case cq_sx:
	case cx_sx:
	case dq_sx:
	case dx_sx:
	case gl_sx:
	case sq_sx:
	case sx_sx:
	case uq_sx:
	case ux_sx:
		return parton::anti_strange;

	case cq_uq:
	case cx_uq:
	case dq_uq:
	case dx_uq:
	case gl_uq:
	case ph_uq:
	case sq_uq:
	case sx_uq:
	case uq_uq:
	case ux_uq:
		return parton::up;

	case cq_ux:
	case dx_ux:
	case sx_ux:
	case uq_ux:
		return parton::anti_up;

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
