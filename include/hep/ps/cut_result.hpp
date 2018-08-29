#ifndef HEP_PS_CUT_RESULT_HPP
#define HEP_PS_CUT_RESULT_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2016-2017  Christopher Schwan
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

#include <utility>

namespace hep
{

/// Class whose instances are used to communicate the result of a phase space
/// cutter back to the class generating the observables, in particular the
/// distributions. The minimum information needed from the cuts needed is
/// whether a given phase space point passes the cuts or not.
template <typename I>
class cut_result_with_info
{
public:
    /// Constructor.
    template <typename NegativeInfo, typename PositiveInfo>
    cut_result_with_info(
        bool neg_cutted,
        bool pos_cutted,
        NegativeInfo&& neg_info,
        PositiveInfo&& pos_info
    )
        : neg_cutted_(neg_cutted)
        , pos_cutted_(pos_cutted)
        , neg_info_(std::forward<NegativeInfo>(neg_info))
        , pos_info_(std::forward<PositiveInfo>(pos_info))
    {
    }

    /// Constructor. Requires the types `N` and `P` to be default constructible.
    cut_result_with_info(bool neg_cutted = true, bool pos_cutted = true)
        : neg_cutted_(neg_cutted)
        , pos_cutted_(pos_cutted)
        , neg_info_()
        , pos_info_()
    {
    }

    /// Type definition for the template parameter `I`.
    using info_t = I;

    /// Returns `true` if a given phase space point has been cutted for negative
    /// rapidity.
    bool neg_cutted() const
    {
        return neg_cutted_;
    }

    /// Returns `true` if a given phase space point has been cutted for positive
    /// rapidity.
    bool pos_cutted() const
    {
        return pos_cutted_;
    }

    /// Returns information about the cuts applied with negative rapidity.
    I const& neg_info() const
    {
        return neg_info_;
    }

    /// Returns information about the cuts applied with positive rapidity.
    I const& pos_info() const
    {
        return pos_info_;
    }

private:
    bool neg_cutted_;
    bool pos_cutted_;
    I neg_info_;
    I pos_info_;
};

/// If this class is used as the template parameter `I` of \ref
/// cut_result_with_info this shows that there is no additional information.
class trivial_cut_info
{
public:
    /// This class has no constructor.
    trivial_cut_info() = delete;
};

/// \cond DOXYGEN_IGNORE
template <>
class cut_result_with_info<trivial_cut_info>
{
public:
    constexpr cut_result_with_info(
        bool neg_cutted = true,
        bool pos_cutted = true
    )
        : neg_cutted_(neg_cutted)
        , pos_cutted_(pos_cutted)
    {
    }

    using info_t = trivial_cut_info;

    constexpr bool neg_cutted() const
    {
        return neg_cutted_;
    }

    constexpr bool pos_cutted() const
    {
        return pos_cutted_;
    }

private:
    bool neg_cutted_;
    bool pos_cutted_;
};
/// \endcond

/// Backward compatible alias for cut results with no additional information.
using cut_result = cut_result_with_info<trivial_cut_info>;

}

#endif
