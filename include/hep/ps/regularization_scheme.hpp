#ifndef HEP_PS_REGULARIZATION_SCHEME_HPP
#define HEP_PS_REGULARIZATION_SCHEME_HPP

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

namespace hep
{

/// Enumeration that denotes the regularization scheme and specifies how the finite parts of the
/// calculation are distributed across the virtual and integrated subtraction terms.
enum class regularization_scheme
{
    /// Dimensional regularization using the Binoth Les Houches Accord (BLHA) convention. This uses
    /// the insertion operator defined in the original Catani-Seymour paper and the \f$
    /// \overline{\mathrm{MS}} \f$ scale replacement \f$ \mu^2 \to \mu^2 \frac{\exp
    /// (\gamma_\text{E})}{4 \pi} \f$.
    dim_reg_blha,

    /// Dimensional regularization using the COLI convention. This differs to \ref
    /// regularization_scheme::dim_reg_blha by replacing the factor \f$ \frac{1}{\Gamma (1-
    /// \epsilon)} \f$ with \f$ \Gamma (1 + \epsilon) \f$ which changes also terms in the constant
    /// part of the Laurent-series.
    dim_reg_coli
};

}

#endif
