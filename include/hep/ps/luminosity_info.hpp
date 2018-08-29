#ifndef HEP_PS_LUMINOSITY_INFO_HPP
#define HEP_PS_LUMINOSITY_INFO_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2016-2018  Christopher Schwan
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

namespace hep
{

/// This class captures the kinematical information that is need to convert the
/// phase-space points, which are calculated by the phase-space generators in
/// the center-of-mass frame of the partons, to the center-of-mass frame of the
/// hadrons (the LAB frame).
template <typename T>
class luminosity_info
{
public:
    /// Constructor.
    luminosity_info(T energy_squared = T())
        : x1_{T(1.0)}
        , x2_{T(1.0)}
        , energy_squared_{energy_squared}
    {
    }

    /// Constructor. See \ref x1, \ref x2, \ref energy_squared, and \ref
    /// rapidity_shift for a description of the arguments.
    luminosity_info(T x1, T x2, T energy_squared, T rapidity_shift)
        : x1_{x1}
        , x2_{x2}
        , energy_squared_{energy_squared}
        , rapidity_shift_{rapidity_shift}
    {
    }

    /// Returns the momentum fraction of the parton coming from the first
    /// hadron.
    T x1() const
    {
        return x1_;
    }

    /// Returns the momentum fraction of the parton coming from the second
    /// hadron.
    T x2() const
    {
        return x2_;
    }

    /// Returns the squared energy of the partons in their center-of-mass frame.
    T energy_squared() const
    {
        return energy_squared_;
    }

    /// Returns the rapidity by which one has to boost the center-of-mass
    /// phase-space point along the beam-axis to get the LAB frame phase-space
    /// point.
    T rapidity_shift() const
    {
        return rapidity_shift_;
    }

private:
    T x1_;
    T x2_;
    T energy_squared_;
    T rapidity_shift_;
};

}

#endif
