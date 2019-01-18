#ifndef HEP_PS_NFPC_HPP
#define HEP_PS_NFPC_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017-2019  Christopher Schwan
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
#include <vector>

namespace hep
{

/// Captures the phase space independent and process specific information to calculate the
/// non-factorizable photonic correction factor in mass regularization.
template <typename T>
struct nfpc_info
{
    /// Mass for each resonance.
    std::vector<T> masses;

    /// For each resonance the vector must contain the indices of the external particles it decays
    /// to.
    std::vector<std::vector<std::size_t>> decays;

    /// The electric charges of all external particles measured in units of the elementary charge.
    /// The charge is the charge of the corresponding particle, irrespectively of whether it is a
    /// particle or an antiparticle.
    std::vector<T> charges;

    /// Signs for all external particles. Incoming antiparticles and outgoing particles are
    /// represented with `-1`, outgoing antiparticles and incoming particles are represented with
    /// `+1`.
    std::vector<T> signs;

    /// Width for each resonance.
    std::vector<T> widths;

    /// The indices of the non-resonant particles that are not decay products of resonances.
    std::vector<std::size_t> non_resonant_particles;

    /// The (unsquared) mass of the photon to mass regularize IR singularities. If the Collier
    /// parameters \f$ \Delta^{(1)}_\text{IR} \f$ and \f$ \mu_\text{IR} \f$ are given, the photon
    /// mass \f$ m_\gamma \f$ is given by the relation
    /// \f[
    ///   m_\gamma^2 = \mu_\text{IR} exp \left( \Delta^{(1)}_\text{IR} \right).
    /// \f]
    T photon_mass;
};

/// Calculates the non-factorizable photonic correction factor in mass regularization for a phase
/// space point in `on_shell_phase_space` and the process specified in `info`. The phase space point
/// must be such that the intermediate resonances given in `info` are exactly on-shell. The vector
/// `resonance_invariants` must contain the virtualities of the resonances for the off-shell phase
/// space point.
template <typename T>
T nfpc(
    std::vector<T> const& on_shell_phase_space,
    std::vector<T> const& resonance_invariants,
    nfpc_info<T> const& info
);

}

#endif
