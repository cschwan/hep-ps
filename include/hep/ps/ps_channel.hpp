#ifndef HEP_PS_PS_CHANNEL_HPP
#define HEP_PS_PS_CHANNEL_HPP

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

#include "hep/ps/ps_decay.hpp"
#include "hep/ps/ps_invariant.hpp"
#include "hep/ps/ps_tchannel.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace hep
{

/// Representation of a phase space channel, a prescription used to generate phase from. A phase
/// space channel is a collection of \ref invariants, \ref tchannels, and \ref decays. Using these
/// three structures phase space is generated according to the following formula:
/// \f[
///     \left( \prod_{i=1}^{n-4} \int_{s_\mathrm{min}^i}^{s_\mathrm{max}^i} \mathrm{d} s_i \right)
///     \left( \prod_{j=1}^{n_\mathrm{t}} \int_{t_\mathrm{min}^j}^{t_\mathrm{max}^j} \mathrm{d} t
///            \int_0^\pi \mathrm{d} \theta_j \right)
///     \left( \prod_{k=1}^{n_\mathrm{s}} \int_0^\pi \mathrm{d} \theta_k
///             \int_0^{2 \pi} \mathrm{d} \varphi_k \right)
/// \f]
/// where \f$ \{ s_i \}_{i=1}^{n-4} \f$ are (timelike-)invariants which are generated in ascending
/// order from their domains \f$ \{ [ s_\mathrm{min}^i, s_\mathrm{max}^i ) \}_{i=1}^{n-4} \f$. The
/// ordering matters, since the domain of the invariant \f$ s_i \f$ may depend on the values of the
/// previous invariants \f$ \{ s_j \}_{j=1}^{i-1} \f$.
class ps_channel
{
public:
    /// Constructor.
    ps_channel(
        std::vector<ps_invariant>&& invariants,
        std::vector<ps_tchannel>&& tchannels,
        std::vector<ps_decay>&& decays
    )
        : invariants_{std::move(invariants)}
        , tchannels_{std::move(tchannels)}
        , decays_{std::move(decays)}
    {
    }


    /// Returns all timelike invariants of this phase space channel.
    std::vector<ps_invariant> const& invariants() const
    {
        return invariants_;
    }

    /// Returns all timelike invariants of this phase space channel.
    std::vector<ps_invariant>& invariants()
    {
        return invariants_;
    }

    /// Returns the spacelike invariants of this phase space channel.
    std::vector<ps_tchannel> const& tchannels() const
    {
        return tchannels_;
    }

    /// Returns the decays of this phase space channel.
    std::vector<ps_decay> const& decays() const
    {
        return decays_;
    }

private:
    std::vector<ps_invariant> invariants_;
    std::vector<ps_tchannel> tchannels_;
    std::vector<ps_decay> decays_;
};

}

#endif
