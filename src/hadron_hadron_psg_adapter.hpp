#ifndef HEP_PS_HH_PHASE_SPACE_GENERATOR_HPP
#define HEP_PS_HH_PHASE_SPACE_GENERATOR_HPP

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

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace hep
{

/// A phase space generator adapter for hadron-hadron collisions. This class
/// takes another phase space generator `G` to generate the momenta of the
/// partons while the class takes care of the additional integration over the
/// momentum fractions.
template <typename G>
class hadron_hadron_psg_adapter
    : public phase_space_generator<typename G::numeric_type>
{
public:
    /// Numeric type used for phase space computations.
    using numeric_type = typename G::numeric_type;

    /// Constructor. The argument `min_energy` is the minimal energy that can be
    /// generated for the partons and `args...` are the arguments for the
    /// constructor of `G`.
    template <typename... Args>
    hadron_hadron_psg_adapter(
        numeric_type min_energy,
        numeric_type cmf_energy,
        Args&&... args
    )
        : psg(std::forward<Args>(args)...)
        , min_energy_(min_energy)
        , cmf_energy_(cmf_energy)
    {
    }

    /// Returns the number of channels of `G`.
    std::size_t channels() const override
    {
        return psg.channels();
    }

    /// Write the densities for each channel of the last generated phase space
    /// point into `densities` and returns an additional jacobian.
    numeric_type densities(std::vector<numeric_type>& densities) override
    {
        return psg.densities(densities) * jacobian_;
    }

    /// Returns the number of dimensions of `G` plus two, which are needed to
    /// generate the momentum fractions.
    std::size_t dimensions() const override
    {
        return psg.dimensions() + 2;
    }

    /// Generates a phase space point.
    void generate(
        std::vector<numeric_type> const& random_numbers,
        std::vector<numeric_type>& momenta,
        std::size_t channel
    ) override {
        using std::log;
        using std::pow;
        using std::sqrt;

        auto const r1 = random_numbers.at(psg.dimensions() + 0);
        auto const r2 = random_numbers.at(psg.dimensions() + 1);

        auto const s = cmf_energy_ * cmf_energy_;
        auto const tau0 = (min_energy_ * min_energy_) / s;
        auto const tau = pow(tau0, r1);
        auto const y = pow(tau, numeric_type(1.0) - r2);
        auto const x1 = y;
        auto const x2 = tau / y;
        auto const shat = tau * s;
        auto const log_tau0 = log(tau0);
        auto const rapidity_shift = log_tau0 * r1 * (r2 - numeric_type(0.5));
        auto const energy = sqrt(shat);

        jacobian_ = r1 * log_tau0 * log_tau0 * tau;
        info_ = luminosity_info<numeric_type>(x1, x2, shat, rapidity_shift);
        psg.generate(random_numbers, momenta, energy, channel);
    }

    /// Returns an instance of \ref luminosity_info that captures the data from
    /// the last generated point.
    luminosity_info<numeric_type> info() const override
    {
        return info_;
    }

    /// Returns `map_dimensions()` of `PhaseSpaceGenerator`.
    std::size_t map_dimensions() const override
    {
        return psg.map_dimensions();
    }

private:
    G psg;
    numeric_type min_energy_;
    numeric_type cmf_energy_;
    numeric_type jacobian_;
    luminosity_info<numeric_type> info_;
};

}

#endif
