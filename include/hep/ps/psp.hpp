#ifndef HEP_PS_PSP_HPP
#define HEP_PS_PSP_HPP

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

#include "hep/ps/psp_type.hpp"
#include "hep/ps/recombined_state.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// Instances of this class represent a single phase space point. Unless otherwise noted the member
/// functions return observables in the laboratory frame.
template <typename T>
class psp
{
public:
    /// Constructor.
    psp(std::vector<T>& cms_psp, std::vector<recombined_state>& states, T rap_shift, psp_type type);

    /// Copy construction is not allowed.
    psp(psp<T> const&) = delete;

    /// Move construction is not allowed.
    psp(psp<T>&&) = delete;

    /// Copy assignment is not allowed.
    psp<T>& operator=(psp<T> const&) = delete;

    /// Move assignment is not allowed.
    psp<T>& operator=(psp<T>&&) = delete;

    /// Calculates the difference of the azimuthal angles of two particles with indices `i` and `j`.
    T abs_phi_diff(std::size_t i, std::size_t j) const;

    /// Returns the cosine of the angle between two particles with indices `i` and `j`.
    T cos_angle(std::size_t i, std::size_t j) const;

    /// Calculates the squared R-distance for two particles with indices `i` and `j`.
    T dist2(std::size_t i, std::size_t j) const;

    /// Returns the invariant mass for the particle with index `i`.
    T m2(std::size_t i) const;

    /// Returns the invariant mass for the particles with indices `i` and `j`.
    T m2(std::size_t i, std::size_t j) const;

    /// Returns the invariant mass for the particles with indices `i`, `j`, and `k`.
    T m2(std::size_t i, std::size_t j, std::size_t k) const;

    /// Returns the invariant mass for four particles with indices `i`, `j`, 'k' and 'l'.
    T m2(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;

    /// Returns the invariant mass for four particles with indices `i`, `j`, 'k', 'l', and `m`.
    T m2(std::size_t i, std::size_t j, std::size_t k, std::size_t l, std::size_t m) const;

    /// Returns the invariant mass for six particles with indices `i`, `j`, 'k', 'l', `m`, and `n`.
    T m2(std::size_t i, std::size_t j, std::size_t k, std::size_t l, std::size_t m, std::size_t n) const;

    /// Returns the transverse mass for the particles with indices `i` and `j`.
    T mt(std::size_t i, std::size_t j) const;

    /// Returns the transverse mass for the particles with indices `i`, `j`, `k`, and `l`.
    T mt(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;

    /// Calculates the squared R-distance for two particles with indices `i` and `j` with
    /// pseudo-rapidity.
    T pdist2(std::size_t i, std::size_t j) const;

    /// Calculates the azimuthal angle for the particle with index `i`.
    T phi(std::size_t i) const;

    /// Calculates the squared transverse momentum of the particle with index `i`.
    T pt2(std::size_t i) const;

    /// Calculates the squared transverse momentum of the particles with indices `i` and `j`.
    T pt2(std::size_t i, std::size_t j) const;

    /// Calculates the squared transverse momentum of the particles with indices `i`, `j`, and `k`.
    T pt2(std::size_t i, std::size_t j, std::size_t k) const;

    /// Calculates the squared transverse momentum of the particles with indices `i`, `j`, `k`, `l`,
    /// `m`, and `n`.
    T pt2(std::size_t i, std::size_t j, std::size_t k, std::size_t l, std::size_t m, std::size_t n) const;

    /// Returns the difference of rapidities of the particles with indices `i` and `j`.
    T rap_diff(std::size_t i, std::size_t j) const;

    /// Returns the difference of rapidities of the particles with indices `i` and `j`.
    T prap_diff(std::size_t i, std::size_t j) const;

    /// Calculates the rapidity of the particle with index `i`.
    T rap(std::size_t i) const;

    /// Calculate the pseudo-rapidity of the particle with index `i`.
    T prap(std::size_t i) const;

    /// Returns the center-of-mass phase space point.
    std::vector<T> const& cms_psp() const;

    /// Returns the rapidity shift.
    T rap_shift() const;

    /// Returns the type of phase space point.
    psp_type type() const;

    /// Removes the object with index `i`.
    void remove(std::size_t i);

    /// Returns the state of the object with index `i`.
    recombined_state state(std::size_t i) const;

    /// Returns all states.
    std::vector<recombined_state> const& states() const;

private:
    std::vector<T>& p_;
    std::vector<recombined_state>& states_;
    T rap_shift_;
    T sign_;
};

}

#endif
