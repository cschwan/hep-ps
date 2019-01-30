#ifndef HEP_PS_INT_DIPOLE_HPP
#define HEP_PS_INT_DIPOLE_HPP

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

#include "hep/ps/dipole_vertex.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/particle_type.hpp"

#include <cassert>
#include <cstddef>

namespace hep
{

///
class int_dipole
{
public:
    ///
    int_dipole(std::size_t emitter, std::size_t spectator, dipole_vertex const& vertex)
        : emitter_{emitter}
        , spectator_{spectator}
        , vertex_{vertex}
    {
        int type = (emitter < 2) | ((spectator < 2) << 1);

        switch (type)
        {
        case 0:
            type_ = insertion_term_type::final_final;
            break;

        case 1:
            type_ = insertion_term_type::initial_final;
            break;

        case 2:
            type_ = insertion_term_type::final_initial;
            break;

        case 3:
            type_ = insertion_term_type::initial_initial;
            break;
        }
    }

    /// Constructor for `born` insertion terms.
    int_dipole(std::size_t initial, dipole_vertex const& vertex)
        : type_{insertion_term_type::born}
        , emitter_{initial}
        , spectator_{0}
        , vertex_{vertex}
    {
    }

    /// Returns the index of the emitter particle.
    std::size_t emitter() const
    {
        assert( type_ != insertion_term_type::born );

        return emitter_;
    }

    /// Returns the index of the spectator particle.
    std::size_t spectator() const
    {
        assert( type_ != insertion_term_type::born );

        return spectator_;
    }

    /// Returns the type of this insertion term.
    insertion_term_type type() const
    {
        return type_;
    }

    /// Returns the index of the initial particle.
    std::size_t initial_particle() const
    {
        switch (type_)
        {
        case insertion_term_type::born:
            return emitter_;

        case insertion_term_type::initial_final:
        case insertion_term_type::initial_initial:
            return emitter();

        case insertion_term_type::final_initial:
            return spectator();

        case insertion_term_type::final_final:
            assert( false );

        default:
            assert( false );
        }
    }

    /// Returns the splitting represented by this integrated dipole.
    dipole_vertex const& vertex() const
    {
        return vertex_;
    }

private:
    insertion_term_type type_;
    std::size_t emitter_;
    std::size_t spectator_;
    dipole_vertex vertex_;
};

/// Comparison operator for two instances of \ref int_dipole `a` and `b`.
bool operator<(int_dipole const& a, int_dipole const& b);

/// Comparison operator for two instances of \ref int_dipole `a` and `b`.
bool operator==(int_dipole const& a, int_dipole const& b);

}

namespace std
{

template <>
struct hash<hep::int_dipole>
{
    size_t operator()(hep::int_dipole const& dipole) const;
};

}

#endif
