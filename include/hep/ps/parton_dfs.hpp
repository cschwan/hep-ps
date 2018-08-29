#ifndef HEP_PS_PARTON_DFS_HPP
#define HEP_PS_PARTON_DFS_HPP

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

#include "hep/ps/parton.hpp"
#include "hep/ps/scales.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace hep
{

/// Class that gives access to parton distribution functions (PDFs); depending
/// on how they are constructed, this class gives access to either the central
/// PDF or the entire PDF set for calculating a central prediction and its
/// uncertainty.
template <typename T>
class parton_dfs
{
public:
    /// Constructor. This creates a single PDF with index `pdf_member` from the
    /// set of PDFs called `name`.
    parton_dfs(std::string const& name, std::size_t pdf_member);

    /// Constructor. This creates all PDFs from the set of PDFs called `name`.
    parton_dfs(std::string const& name);

    /// Move constructor.
    parton_dfs(parton_dfs<T>&& pdf);

    /// Destructor.
    ~parton_dfs();

    /// Evaluates the strong coupling for the renormalization scales in `scales`
    /// and writes the results into `alphas`
    void eval_alphas(
        std::vector<hep::scales<T>> const& scales,
        std::vector<T>& alphas
    );

    /// This function returns the number of PDFs represented with by this
    /// object.
    std::size_t count() const;

    /// Evaluates the central PDFs for all partons for the momentum fraction `x`
    /// at the given `scales`, writes the results into `scale_pdfs`, and also
    /// evaluates, for the central scale, the uncertainty PDFs and writes the
    /// results into `uncertainty_pdfs`.
    void eval(
        T x,
        std::vector<scales<T>> const& scales,
        std::vector<parton_array<T>>& scale_pdfs,
        std::vector<parton_array<T>>& uncertainty_pdfs
    );

    /// Registers all the partons contained in `set`. For performance reasons
    /// the values for the PDFs of unregistered partons may be zero.
    void register_partons(parton_set set);

private:
    class impl;
    std::unique_ptr<impl> pimpl;
};

}

#endif
