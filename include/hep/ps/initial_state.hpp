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
	bq_cq, /*!< \f$ Q=1/3 \f$ bottom charm initial state. */
	bq_dx, /*!< \f$ Q=0 \f$ bottom anti-down initial state. */
	bq_sx, /*!< \f$ Q=0 \f$ bottom anti-strange initial state. */
	bq_uq, /*!< \f$ Q=1/3 \f$ bottom up initial state. */
	bq_gl, /*!< \f$ Q=-1/3 \f$ bottom gluon initial state. */
	bx_cq, /*!< \f$ Q=1 \f$ anti-bottom charm initial state. */
	bx_dx, /*!< \f$ Q=2/3 \f$ anti-bottom anti-down initial state. */
	bx_sx, /*!< \f$ Q=2/3 \f$ anti-bottom anti-strange initial state. */
	bx_uq, /*!< \f$ Q=1 \f$ anti-bottom up initial state. */
	bx_gl, /*!< \f$ Q=1/3 \f$ anti-bottom gluon initial state. */
	cq_bq, /*!< \f$ Q=1/3 \f$ charm bottom initial state. */
	cq_bx, /*!< \f$ Q=1 \f$ charm anti-bottom initial state. */
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
	cx_gl, /*!< \f$ Q=-2/3 \f$ anti-charm gluon initial state. */
	cx_ph, /*!< \f$ Q=-2/3 \f$ anti-charm photon initial state. */
	cx_sx, /*!< \f$ Q=-1/3 \f$ anti-charm anti-strange initial state. */
	cx_uq, /*!< \f$ Q=0 \f$ anti-charm up initial state. */
	dq_cq, /*!< \f$ Q=1/3 \f$ down charm initial state. */
	dq_dx, /*!< \f$ Q=0 \f$ down anti-down initial state. */
	dq_gl, /*!< \f$ Q=-1/3 \f$ down gluon initial state. */
	dq_ph, /*!< \f$ Q=-1/3 \f$ down photon initial state. */
	dq_sx, /*!< \f$ Q=0 \f$ down anti-strange initial state. */
	dq_uq, /*!< \f$ Q=1/3 \f$ down up initial state. */
	dx_bq, /*!< \f$ Q=0 \f$ anti-down bottom initial state. */
	dx_bx, /*!< \f$ Q=2/3 \f$ anti-down anti-bottom initial state. */
	dx_cq, /*!< \f$ Q=1 \f$ anti-down charm initial state. */
	dx_cx, /*!< \f$ Q=-1/3 \f$ anti-down anti-charm initial state. */
	dx_dq, /*!< \f$ Q=0 \f$ anti-down down initial state. */
	dx_dx, /*!< \f$ Q=2/3 \f$ anti-down anti-down initial state. */
	dx_gl, /*!< \f$ Q=1/3 \f$ anti-down gluon initial state. */
	dx_ph, /*!< \f$ Q=1/3 \f$ anti-down photon initial state. */
	dx_sq, /*!< \f$ Q=0 \f$ anti-down strange initial state. */
	dx_sx, /*!< \f$ Q=2/3 \f$ anti-strange anti-down initial state. */
	dx_uq, /*!< \f$ Q=1 \f$ anti-down up initial state. */
	dx_ux, /*!< \f$ Q=-1/3 \f$ anti-down anti-up initial state. */
	gl_bq, /*!< \f$ Q=-1/3 \f$ gluon bottom initial state */
	gl_bx, /*!< \f$ Q=1/3 \f$ gluon anti-bottom initial state */
	gl_cq, /*!< \f$ Q=2/3 \f$ gluon charm initial state. */
	gl_cx, /*!< \f$ Q=-2/3 \f$ gluon anti-charm initial state. */
	gl_dq, /*!< \f$ Q=1/3 \f$ gluon down initial state. */
	gl_dx, /*!< \f$ Q=1/3 \f$ gluon anti-down initial state. */
	gl_gl, /*!< \f$ Q=0 \f$ gluon gluon initial state. */
	gl_ph, /*!< \f$ Q=0 \f$ gluon photon initial state. */
	gl_sq, /*!< \f$ Q=-1/3 \f$ gluon strange initial state. */
	gl_sx, /*!< \f$ Q=1/3 \f$ gluon anti-strange initial state. */
	gl_uq, /*!< \f$ Q=2/3 \f$ gluon up initial state. */
	gl_ux, /*!< \f$ Q=-2/3 \f$ gluon anti-up initial state. */
	ph_cq, /*!< \f$ Q=2/3 \f$ photon charm initial state. */
	ph_dx, /*!< \f$ Q=1/3 \f$ photon anti-down initial state. */
	ph_gl, /*!< \f$ Q=0 \f$ photon gluon initial state. */
	ph_ph, /*!< \f$ Q=0 \f$ photon photon initial state. */
	ph_sx, /*!< \f$ Q=1/3 \f$ photon anti-strange initial state. */
	ph_uq, /*!< \f$ Q=2/3 \f$ photon up initial state. */
	sq_cq, /*!< \f$ Q=1/3 \f$ strange charm initial state. */
	sq_dx, /*!< \f$ Q=0 \f$ strange anti-down initial state. */
	sq_gl, /*!< \f$ Q=-1/3 \f$ strange gluon initial state. */
	sq_ph, /*!< \f$ Q=-1/3 \f$ strange photon initial state. */
	sq_sx, /*!< \f$ Q=0 \f$ strange anti-strange initial state. */
	sq_uq, /*!< \f$ Q=1/3 \f$ strange up initial state. */
	sx_bq, /*!< \f$ Q=0 \f$ anti-strange bottom initial state. */
	sx_bx, /*!< \f$ Q=2/3 \f$ anti-strange anti-bottom initial state. */
	sx_cq, /*!< \f$ Q=1 \f$ anti-strange charm initial state. */
	sx_cx, /*!< \f$ Q=-1/3 \f$ anti-strange anti-charm initial state. */
	sx_dq, /*!< \f$ Q=0 \f$ anti-strange down initial state. */
	sx_dx, /*!< \f$ Q=2/3 \f$ anti-strange anti-down initial state. */
	sx_gl, /*!< \f$ Q=1/3 \f$ anti-strange gluon initial state. */
	sx_ph, /*!< \f$ Q=1/3 \f$ anti-strange photon initial state. */
	sx_sq, /*!< \f$ Q=0 \f$ anti-strange strange initial state. */
	sx_sx, /*!< \f$ Q=2/3 \f$ anti-strange anti-strange initial state. */
	sx_uq, /*!< \f$ Q=1 \f$ anti-strange up initial state. */
	sx_ux, /*!< \f$ Q=-1/3 \f$ anti-strange anti-up initial state. */
	uq_bq, /*!< \f$ Q=1/3 \f$ up bottom initial state. */
	uq_bx, /*!< \f$ Q=1 \f$ up anti-bottom initial state. */
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
	ux_gl, /*!< \f$ Q=-2/3 \f$ anti-up gluon initial state. */
	ux_ph, /*!< \f$ Q=-2/3 \f$ anti-up photon initial state. */
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
	case bq_cq:
	case bq_dx:
	case bq_gl:
	case bq_sx:
	case bq_uq:
		return parton::bottom;

	case bx_cq:
	case bx_dx:
	case bx_gl:
	case bx_sx:
	case bx_uq:
		return parton::anti_bottom;

	case cq_bq:
	case cq_bx:
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
	case cx_gl:
	case cx_ph:
	case cx_sx:
	case cx_uq:
		return parton::anti_charm;

	case dq_cq:
	case dq_dx:
	case dq_gl:
	case dq_ph:
	case dq_sx:
	case dq_uq:
		return parton::down;

	case dx_bq:
	case dx_bx:
	case dx_cq:
	case dx_cx:
	case dx_dq:
	case dx_dx:
	case dx_gl:
	case dx_ph:
	case dx_sq:
	case dx_sx:
	case dx_uq:
	case dx_ux:
		return parton::anti_down;

	case gl_bq:
	case gl_bx:
	case gl_cq:
	case gl_cx:
	case gl_dq:
	case gl_dx:
	case gl_gl:
	case gl_ph:
	case gl_sq:
	case gl_sx:
	case gl_uq:
	case gl_ux:
		return parton::gluon;

	case ph_cq:
	case ph_dx:
	case ph_gl:
	case ph_ph:
	case ph_sx:
	case ph_uq:
		return parton::photon;

	case sq_cq:
	case sq_dx:
	case sq_gl:
	case sq_ph:
	case sq_sx:
	case sq_uq:
		return parton::strange;

	case sx_bq:
	case sx_bx:
	case sx_cq:
	case sx_cx:
	case sx_dq:
	case sx_dx:
	case sx_gl:
	case sx_ph:
	case sx_sq:
	case sx_sx:
	case sx_uq:
	case sx_ux:
		return parton::anti_strange;

	case uq_bq:
	case uq_bx:
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
	case ux_gl:
	case ux_ph:
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
	case gl_bq:
	case cq_bq:
	case dx_bq:
	case sx_bq:
	case uq_bq:
		return parton::bottom;

	case gl_bx:
	case cq_bx:
	case dx_bx:
	case sx_bx:
	case uq_bx:
		return parton::anti_bottom;

	case bq_cq:
	case bx_cq:
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
	case gl_cx:
	case sx_cx:
	case uq_cx:
		return parton::anti_charm;

	case cq_dq:
	case dx_dq:
	case gl_dq:
	case sx_dq:
	case uq_dq:
		return parton::down;

	case bq_dx:
	case bx_dx:
	case cq_dx:
	case cx_dx:
	case dq_dx:
	case dx_dx:
	case gl_dx:
	case ph_dx:
	case sq_dx:
	case sx_dx:
	case uq_dx:
	case ux_dx:
		return parton::anti_down;

	case bq_gl:
	case bx_gl:
	case cq_gl:
	case cx_gl:
	case dq_gl:
	case dx_gl:
	case gl_gl:
	case ph_gl:
	case sq_gl:
	case sx_gl:
	case uq_gl:
	case ux_gl:
		return parton::gluon;

	case cq_ph:
	case cx_ph:
	case dq_ph:
	case dx_ph:
	case gl_ph:
	case ph_ph:
	case sq_ph:
	case sx_ph:
	case uq_ph:
	case ux_ph:
		return parton::photon;

	case cq_sq:
	case dx_sq:
	case gl_sq:
	case sx_sq:
	case uq_sq:
		return parton::strange;

	case bq_sx:
	case bx_sx:
	case cq_sx:
	case cx_sx:
	case dq_sx:
	case dx_sx:
	case gl_sx:
	case ph_sx:
	case sq_sx:
	case sx_sx:
	case uq_sx:
	case ux_sx:
		return parton::anti_strange;

	case bq_uq:
	case bx_uq:
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
	case gl_ux:
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

/// Returns an `initial_state` constructed from two \ref parton instances.
initial_state partons_to_initial_state(parton one, parton two);

/// Returns the value of the Casimir operator for the particle in the \ref
/// initial_state `state` with the given `index`.
template <typename T>
T casimir_operator(initial_state state, std::size_t index);

}

#endif
