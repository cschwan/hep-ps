#ifndef HEP_PS_CONVOLUTE_HPP
#define HEP_PS_CONVOLUTE_HPP

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

#include "hep/ps/initial_state.hpp"
#include "hep/ps/parton.hpp"

#include "nonstd/span.hpp"

#include <cassert>
#include <vector>

namespace hep
{

template <typename T>
inline T convolute(
    parton_array<T> const& pdf_x1,
    parton_array<T> const& pdf_x2,
    initial_state_map<T> const& matrix_elements,
    initial_state_set set
) {
    T result = T();

    for (auto const& keyvals : matrix_elements)
    {
        auto const state = keyvals.first;

        if (!set.includes(state))
        {
            continue;
        }

        auto const one = state_parton_one(state);
        auto const two = state_parton_two(state);
        auto const sym = (one == two) ? T(0.5) : T(1.0);
        auto const me = keyvals.second;

        result += sym * pdf_x1[one] * pdf_x2[two] * me;
    }

    return result;
}

template <typename T>
inline void convolute_mes_with_pdfs(
    std::vector<T>& scale_results,
    std::vector<T>& pdf_results,
    nonstd::span<parton_array<T> const> scale_uncertainty_pdfs_one,
    nonstd::span<parton_array<T> const> scale_uncertainty_pdfs_two,
    nonstd::span<parton_array<T> const> pdf_uncertainty_pdfs_one,
    nonstd::span<parton_array<T> const> pdf_uncertainty_pdfs_two,
    nonstd::span<initial_state_map<T> const> matrix_elements,
    initial_state_set set,
    nonstd::span<T const> alphas_factors,
    T global_factor
) {
    std::size_t const scales = scale_uncertainty_pdfs_one.size();
    std::size_t const pdfs = pdf_uncertainty_pdfs_one.size();

    assert( size(scale_uncertainty_pdfs_one) == scales );
    assert( size(scale_uncertainty_pdfs_two) == scales );
    assert( size(matrix_elements) == scales );
    assert( size(alphas_factors) == scales );
    assert( size(pdf_uncertainty_pdfs_one) == pdfs );
    assert( size(pdf_uncertainty_pdfs_two) == pdfs );

    scale_results.clear();

    for (std::size_t i = 0; i != scales; ++i)
    {
        scale_results.push_back(global_factor * alphas_factors[i] * convolute(
            scale_uncertainty_pdfs_one[i],
            scale_uncertainty_pdfs_two[i],
            matrix_elements[i],
            set
        ));
    }

    pdf_results.clear();

    for (std::size_t i = 0; i != pdfs; ++i)
    {
        pdf_results.push_back(global_factor * alphas_factors[0] * convolute(
            pdf_uncertainty_pdfs_one[i],
            pdf_uncertainty_pdfs_two[i],
            matrix_elements[0],
            set
        ));
    }
}

}

#endif
