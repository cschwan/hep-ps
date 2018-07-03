#ifndef HEP_PS_PARTON_HPP
#define HEP_PS_PARTON_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017-2018  Christopher Schwan
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

namespace hep
{

HEP_ENUM(parton,
	anti_up,      /*!< An anti-up quark. */
	anti_down,    /*!< An anti-down quark. */
	anti_charm,   /*!< An anti-charm quark. */
	anti_strange, /*!< An anti-strange quark. */
	anti_bottom,  /*!< An anti-bottom quark. */
	gluon,        /*!< A gluon. */
	up,           /*!< An up quark. */
	down,         /*!< A down quark. */
	charm,        /*!< A charm quark. */
	strange,      /*!< A strange quark. */
	bottom,       /*!< A bottom quark. */
	photon        /*!< A photon. */
);

HEP_ENUM_ARRAY(parton);

HEP_ENUM_SET(parton);

}

#endif
