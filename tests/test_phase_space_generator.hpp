#ifndef TEST_PHASE_SPACE_GENERATOR_HPP
#define TEST_PHASE_SPACE_GENERATOR_HPP

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

#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/phase_space_generator.hpp"

#include "catch2/catch.hpp"

#include <cstddef>
#include <vector>

template <typename T>
class test_phase_space_generator : public hep::phase_space_generator<T>
{
public:
    test_phase_space_generator(
        std::size_t final_states,
        std::size_t extra_random_numbers = 0
    )
        : final_states_{final_states}
        , extra_random_numbers_{extra_random_numbers}
    {
    }

    std::size_t channels() const override
    {
        return 1;
    }

    T densities(std::vector<T>& densities) override
    {
        densities.front() = T(1.0);

        return T(1.0);
    }

    std::size_t dimensions() const override
    {
        // we don't need random numbers in this generator
        return 0;
    }

    void generate(
        std::vector<T> const& random_numbers,
        std::vector<T>& momenta,
        std::size_t channel
    ) override {
        REQUIRE( random_numbers.size() == 0 );
        REQUIRE( momenta.size() == map_dimensions() );
        REQUIRE( channel == 0 );

        // fill with zeros
        momenta.assign(map_dimensions(), T{});
    }

    hep::luminosity_info<T> info() const override
    {
        // energy^2 is `0.5` to cancel the factor `2` in the integrands
        return { T{0.5}, T{0.5}, T{0.5}, T{} };
    }

    std::size_t map_dimensions() const override
    {
        return 4 * (final_states_ + 2) + extra_random_numbers_;
    }

private:
    std::size_t final_states_;
    std::size_t extra_random_numbers_;
};

#endif
