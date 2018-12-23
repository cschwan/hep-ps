#ifndef HEP_PS_PS_FUNCTIONS_HPP
#define HEP_PS_PS_FUNCTIONS_HPP

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

#include <array>

namespace hep
{

/// Boosts the vector `p` in the direction given by `q`, whose invariant must be `m`. If `inverse`
/// is `true` the inverse boost is applied, otherwise the (uninverted) boost, which t ransforms `p`
/// to the rest frame of `q`.
template <typename T>
void boost(T m, std::array<T, 4>& p, std::array<T, 4> const& q, bool inverse);

/// Calculates the Källén function, also known as triangle function, that arises in three-point
/// one-loop functions or particle decays.
template <typename T>
T kaellen(T x, T y, T z);

/// Calculates the square root of the Källén function.
template <typename T>
T sqrt_kaellen(T x, T y, T z);

/// Rotates the vector `p` first around the 3-axis by `theta` and then around the 2-axis by `phi`.
template <typename T>
void rotate(std::array<T, 4>& p, T phi, T cos_theta);

/// Captures the result of the function \ref calculate_space_like_invariant_bounds.
template <typename T>
struct space_like_invariant_bounds
{
    /// Minimum value of the space-like invariant.
    T min;

    /// Maximum value of the space-like invariant.
    T max;

    /// Square-root of the Källén function \f$ \sqrt{\lambda_s} \f$ of `s` with the remaining
    /// light/time-like invariants: \f$ \lambda_s = \lambda (s, s_1, s_2) \f$.
    T lambdas;

    /// Square-root of the Källén function \f$ \sqrt{\lambda_t} \f$ of `s` with the space-like
    /// invariants: \f$ \lambda_t = \lambda (s, t_1, t_2) \f$.
    T lambdat;
};

/// Calculates the bounds on a space-like invariant determined by the three light/time-like
/// invariants `s`, `s1`, and `s2` and the two space-like invariants `t1` and `t2`.
template <typename T>
space_like_invariant_bounds<T> calculate_space_like_invariant_bounds(T s, T s1, T s2, T t1, T t2);

}

#endif
