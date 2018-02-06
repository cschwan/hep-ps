#ifndef HEP_PS_PARTON_TYPE_HPP
#define HEP_PS_PARTON_TYPE_HPP

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

#include "hep/ps/enum.hpp"
#include "hep/ps/parton.hpp"

namespace hep
{

// TODO: rename `gluon_` to `gluon` once the `enum` is an `enum class`

HEP_ENUM(parton_type,
	anti_quark, /*!< An anti-quark. */
	gluon_,     /*!< A gluon. */
	quark,      /*!< A quark. */
	photon_     /*!< A photon. */
);

HEP_ENUM_ARRAY(parton_type);

/// Returns the \ref parton_type for a \ref parton `p`.
constexpr parton_type parton_type_of(parton p)
{
	switch (p)
	{
	case anti_up:
	case anti_down:
	case anti_charm:
	case anti_strange:
		return parton_type::anti_quark;

	case gluon:
		return parton_type::gluon_;

	case up:
	case down:
	case charm:
	case strange:
		return parton_type::quark;

	case photon:
		return parton_type::photon_;
	}
}

}

#endif
