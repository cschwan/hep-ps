#ifndef HEP_PS_DIPOLE_SPLIT_HPP
#define HEP_PS_DIPOLE_SPLIT_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2019  Christopher Schwan
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

#include "hep/ps/correction_type.hpp"

namespace hep
{

/// Enumeration listing the different splittings that give rise to infrared singularities which are
/// cancelled with the subtraction formulae. The enumeration lists all values which are important
/// for the integrated dipoles, which, in the case of an emitted photon, are sensitive to the charge
/// of the emitted particle.
enum class dipole_split
{
    /// (Anti-)quark splitting into (anti-)quark and a gluon
    quark_to_quark_gluon,

    ///
    gluon_to_gluon_gluon,

    ///
    gluon_to_quark_quark,

    ///
    up_to_up_photon,

    ///
    down_to_down_photon,

    ///
    charm_to_charm_photon,

    ///
    strange_to_strange_photon,

    ///
    photon_to_up_up,

    ///
    photon_to_down_down,

    ///
    photon_to_charm_charm,

    ///
    photon_to_strange_strange
};

///
constexpr correction_type correction_type_of(dipole_split splitting)
{
    switch (splitting)
    {
    case dipole_split::quark_to_quark_gluon:
    case dipole_split::gluon_to_gluon_gluon:
    case dipole_split::gluon_to_quark_quark:
        return correction_type::qcd;

    case dipole_split::up_to_up_photon:
    case dipole_split::down_to_down_photon:
    case dipole_split::charm_to_charm_photon:
    case dipole_split::strange_to_strange_photon:
    case dipole_split::photon_to_up_up:
    case dipole_split::photon_to_down_down:
    case dipole_split::photon_to_charm_charm:
    case dipole_split::photon_to_strange_strange:
        return correction_type::ew;

#if __GNUC__ > 5
    default:
        throw 0;
#endif
    }
}

}

#endif
