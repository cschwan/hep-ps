#ifndef HEP_PS_PARTON_DF_UNC_HPP
#define HEP_PS_PARTON_DF_UNC_HPP

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
#include <memory>
#include <string>
#include <vector>

namespace hep
{

/// Structure that captures the symmetric and asymmetric PDF uncertainties.
template <typename T>
struct pdf_unc
{
    /// Symmetric PDF uncertainties.
    T sym;

    /// Asymmetric PDF uncertainties (lower value).
    T neg;

    /// Asymmetric PDF uncertainties (upper value).
    T pos;
};

/// Class that calculates the PDF uncertainty for a given PDF set.
template <typename T>
class parton_df_unc
{
public:
    /// Constructor.
    parton_df_unc(std::string const& name, T cl = T());

    /// Move constructor.
    parton_df_unc(parton_df_unc&& pdf_unc);

    /// Destructor.
    ~parton_df_unc();

    /// Calculates the PDF uncertainty for the given `values`. The first element
    /// of `values` must be the central prediction, all other values the ones
    /// from the Error PDF set.
    pdf_unc<T> uncertainty(std::vector<T> const& values) const;

private:
    class impl;
    std::unique_ptr<impl> pimpl;
};

}

#endif
