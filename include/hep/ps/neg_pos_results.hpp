#ifndef HEP_PS_NEG_POS_RESULTS_HPP
#define HEP_PS_NEG_POS_RESULTS_HPP

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

namespace hep
{

/**
 * Instances of this class capture the result of (parts of) perturbative
 * calculations. The result is subdivided into a positive and negative rapidity
 * part, which calculated in different reference frames, defined by the
 * transformation rapidities into the LAB frame:
 * \f{align}{
 *     y_\text{neg} &= - y_\text{CMS} + y_\text{PS} \\
 *     y_\text{pos} &= + y_\text{CMS} + y_\text{PS}
 * \f}
 * where \f$ y_\text{CMS} \f$ is the rapidity of any object in the partonic
 * center-of-mass frame and \f$ y_\text{PS} \f$ is the rapdity shift produced
 * by the phase space generator. The results \f$ y_\text{neg}, y_\text{pos} \f$
 * are then the rapidities in the LAB frame.
 */
template <typename T>
struct neg_pos_results
{
    /// Constructor.
    neg_pos_results(T neg_result = T(), T pos_result = T())
        : neg{neg_result}
        , pos{pos_result}
    {
    }

    /// Result for the negative rapidity.
    T neg;

    /// Result for the positive rapidity.
    T pos;
};

/// Addition operator. Adds `other` to `result` and returns `result`.
template <typename T>
inline neg_pos_results<T>& operator+=(
    neg_pos_results<T>& result,
    neg_pos_results<T> const& other
) {
    result.neg += other.neg;
    result.pos += other.pos;

    return result;
}

}

#endif
