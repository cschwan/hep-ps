#ifndef HEP_PS_SPIN_CORRELATION_MATRIX_HPP
#define HEP_PS_SPIN_CORRELATION_MATRIX_HPP

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

/// Class representing a spin-correlation matrix. The matrix is of the form
/// \f[
///     -g^{\mu \nu} A + \frac{p^\mu p^\nu}{p^2} B
/// \f]
/// and instances of this class represent the factors \f$ A \f$, \f$ B \f$, and
/// the vector \f$ p \f$.
template <typename T>
struct spin_correlation_matrix
{
    /// The factor \f$ A \f$.
    T a;

    /// The factor \f$ B \f$.
    T b;

    /// The vector \f$ p \f$.
    std::array<T, 4> p;
};

}

#endif
