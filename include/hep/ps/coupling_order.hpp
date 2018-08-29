#ifndef HEP_PS_COUPLING_ORDER_HPP
#define HEP_PS_COUPLING_ORDER_HPP

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

#include <cstddef>

namespace hep
{

/// Class to identify different coupling orders in the electromagnetic and
/// strong couplings: \f$ \alpha^n \alpha_\mathrm{s}^m \f$. This class saves the
/// powers \f$ n \f$ and \f$ m \f$.
class coupling_order
{
public:
    /// Constructor.
    coupling_order(std::size_t alpha_power, std::size_t alphas_power)
        : alpha_power_{alpha_power}
        , alphas_power_{alphas_power}
    {
    }

    /// Returns the power of \f$ \alpha \f$.
    std::size_t alpha_power() const
    {
        return alpha_power_;
    }

    /// Returns the power of \f$ \alpha_\mathrm{s} \f$.
    std::size_t alphas_power() const
    {
        return alphas_power_;
    }

private:
    std::size_t alpha_power_;
    std::size_t alphas_power_;
};

}

#endif
