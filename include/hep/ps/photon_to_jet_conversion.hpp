#ifndef HEP_PS_PHOTON_TO_JET_CONVERSION_HPP
#define HEP_PS_PHOTON_TO_JET_CONVERSION_HPP

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

#include <cstddef>

namespace hep
{

///
template <typename T>
class photon_to_jet_conversion
{
public:
    /// Default constructor. Constructs a conversion function where \ref active always returns
    /// `false`.
    photon_to_jet_conversion();

    ///
    photon_to_jet_conversion(
        T alpha,
        T delta_alpha_hadr,
        T conversion_scale,
        T nc,
        std::size_t nf
    );

    ///
    bool active() const;

    ///
    T eval(T regularization_scale, T non_perturbative_factor) const;

private:
    T alpha_;
    T delta_alpha_hadr_;
    T conversion_scale_;
    T nc_;
    T nf_;
};

}

#endif
