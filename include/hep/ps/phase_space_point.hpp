#ifndef HEP_PS_PHASE_SPACE_VIEW_HPP
#define HEP_PS_PHASE_SPACE_VIEW_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017-2018  Christopher Schwan
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

/// Class that represents a single phase space point in the center-of-mass frame
/// (CMF) of the (partonic) collision. The member functions of this class enable
/// one to calculate kinematic observables in the LAB frame.
template <typename T>
class phase_space_point
{
public:
    /// Constructor. The phase space point in `phase_space` must be given in the
    /// CM frame. The parameter `rapidity_shift` enables this class to boost
    /// into the LAB frame.
    phase_space_point(
        std::vector<T> const& phase_space,
        T rapidity_shift = T()
    );

    /// Calculates the difference of the azimuthal angles of two particles with
    /// indices `i` and `j`.
    T abs_phi_diff(std::size_t i, std::size_t j) const;

    /// Returns the cosine of the angle between two particles with indices `i`
    /// and `j`.
    T cos_angle_neg(std::size_t i, std::size_t j) const;

    /// Returns the cosine of the angle between two particles with indices `i`
    /// and `j`.
    T cos_angle_pos(std::size_t i, std::size_t j) const;

    /// Calculates the squared R-distance for two particles with indices `i` and
    /// `j`.
    T dist2(std::size_t i, std::size_t j) const;

    /// Returns the invariant mass for the particle with index `i`.
    T m2(std::size_t i) const;

    /// Returns the invariant mass for the particles with indices `i` and `j`.
    T m2(std::size_t i, std::size_t j) const;
    
     /// Returns the invariant mass for 4 particles with indices `i`, `j`, 'k' and 'l'.
    T m2(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;

    /// Returns the transverse mass for the particles with indices `i` and `j`.
    T mt(std::size_t i, std::size_t j) const;

    /// Returns the transverse mass for the particles with indices `i`, `j`,
    /// `k`, and `l`.
    T mt(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const;

    /// Calculates the squared R-distance for two particles with indices `i` and
    /// `j` with pseudo-rapidity.
    T pdist2(std::size_t i, std::size_t j) const;

    /// Calculates the azimuthal angle for the particle with index `i`.
    T phi(std::size_t i) const;

    /// Calculates the squared transverse momentum of the particle with index
    /// `i`.
    T pt2(std::size_t i) const;

    /// Calculates the squared transverse momentum of the particles with indices
    /// `i` and `j`.
    T pt2(std::size_t i, std::size_t j) const;

    /// Returns the difference of rapidities of the particles with indices `i`
    /// and `j`.
    T rap_diff(std::size_t i, std::size_t j) const;

    /// Returns the difference of rapidities of the particles with indices `i`
    /// and `j`.
    T prap_diff(std::size_t i, std::size_t j) const;

    /// Calculates the rapidity of the particle with index `i`.
    T rap_neg(std::size_t i) const;

    /// Calculates the rapidity of the particle with index `i`.
    T rap_pos(std::size_t i) const;

    // Calculate the pseudo-rapidity of the particle with index `i`.
    T prap_neg(std::size_t i) const;

    // Calculate the pseudo-rapidity of the particle with index `i`.
    T prap_pos(std::size_t i) const;

private:
    T const* p_;
    T rapidity_shift_;
};

}

#endif
