#ifndef HEP_PS_PS_DECAY_HPP
#define HEP_PS_PS_DECAY_HPP

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

class ps_decay
{
public:
    ps_decay(std::size_t in_mom_id, std::size_t out1_mom_id, std::size_t out2_mom_id)
        : in_mom_id_{in_mom_id}
        , out1_mom_id_{out1_mom_id}
        , out2_mom_id_{out2_mom_id}
    {
    }

    std::size_t in_mom_id() const
    {
        return in_mom_id_;
    }

    std::size_t out1_mom_id() const
    {
        return out1_mom_id_;
    }

    std::size_t out2_mom_id() const
    {
        return out2_mom_id_;
    }

private:
    std::size_t in_mom_id_;
    std::size_t out1_mom_id_;
    std::size_t out2_mom_id_;
};

}

#endif
