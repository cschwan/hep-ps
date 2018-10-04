#ifndef HEP_PS_STATIC_SCALE_FUNCTION_HPP
#define HEP_PS_STATIC_SCALE_FUNCTION_HPP

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

#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scale_variation.hpp"
#include "hep/ps/scales.hpp"

#include <vector>

namespace hep
{

/// A generic static scale function.
template <typename T>
class static_scale_function
{
public:
    /// Constructor.
    static_scale_function(T scale, scale_variation variation);

    /// Returns the static scale(s).
    void operator()(
        std::vector<T> const& phase_space,
        std::vector<hep::scales<T>>& scales,
        std::vector<recombined_state> const& states
    );

    /// Returns `false`, since this scale function is independent of the phase space point (static).
    bool dynamic() const;

private:
    T scale_;
    scale_variation variation_;
};

}

#endif
