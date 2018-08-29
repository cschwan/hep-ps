#ifndef HEP_PS_SCALE_VARIATION_HPP
#define HEP_PS_SCALE_VARIATION_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2018  Christopher Schwan
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

namespace hep
{

/// Enumeration describing common scale variations.
enum class scale_variation
{
    /// No scale variation.
    single_scale,

    /// A three point variation setting the scales \f$ (\mu_\text{F}, \mu,
    /// \mu_\text{R}) \f$ to the three points \f$ \mu_\text{F} = \mu_\text{R} =
    /// \mu \f$, \f$ \mu_\text{F} = \mu_\text{R} = \mu / 2 \f$, and \f$
    /// \mu_\text{F} = \mu_\text{R} = 2 \mu \f$ while leaving the regularization
    /// scale \f$ \mu \f$ constant.
    three_point,

    /// A seven point variation.
    seven_point
};

}

#endif
