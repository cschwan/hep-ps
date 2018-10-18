#ifndef HEP_PS_CONVOLUTE_HPP
#define HEP_PS_CONVOLUTE_HPP

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

#include "hep/ps/cut_result.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/psp_type.hpp"

#include <cassert>

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

template <typename T, typename I>
inline neg_pos_results<T> convolute(
    parton_array<T> const& pdfx1,
    parton_array<T> const& pdfx2,
    initial_state_map<T> const& matrix_elements,
    initial_state_set set,
    T factor,
    cut_result_with_info<I> const& cut
) {
    T neg{};
    T pos{};

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

        if (!cut.pos_cutted())
        {
            pos += sym * pdfx1[one] * pdfx2[two] * me;
        }

        if (!cut.neg_cutted())
        {
            neg += sym * pdfx1[two] * pdfx2[one] * me;
        }
    }

    return { factor * neg , factor * pos };
}

template <typename T, typename I>
inline neg_pos_results<T> convolute(
    parton_array<T> const& pdfx1_neg,
    parton_array<T> const& pdfx2_neg,
    parton_array<T> const& pdfx1_pos,
    parton_array<T> const& pdfx2_pos,
    initial_state_map<T> const& matrix_elements,
    initial_state_set set,
    T factor,
    cut_result_with_info<I> const& cut
) {
    T neg{};
    T pos{};

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

        if (!cut.pos_cutted())
        {
            pos += sym * pdfx1_pos[one] * pdfx2_pos[two] * me;
        }

        if (!cut.neg_cutted())
        {
            neg += sym * pdfx1_neg[one] * pdfx2_neg[two] * me;
        }
    }

    return { factor * neg , factor * pos };
}

template <typename T>
inline void convolute_mes_with_pdfs(
    psp_type type,
    std::vector<T>& scale_results,
    std::vector<T>& pdf_results,
    std::vector<parton_array<T>> const& scale_uncertainty_pdfs_one,
    std::vector<parton_array<T>> const& scale_uncertainty_pdfs_two,
    std::vector<parton_array<T>> const& pdf_uncertainty_pdfs_one,
    std::vector<parton_array<T>> const& pdf_uncertainty_pdfs_two,
    std::vector<initial_state_map<T>> const& matrix_elements,
    initial_state_set set,
    std::vector<T> const& alphas_factors,
    T global_factor
) {
    std::size_t const scale_index = (type == psp_type::pos_rap) ? 0 :
        scale_uncertainty_pdfs_one.size() / 2;
    std::size_t const scales = scale_uncertainty_pdfs_one.size() / 2;
    std::size_t const pdf_index = (type == psp_type::pos_rap) ? 0 :
        pdf_uncertainty_pdfs_one.size() / 2;
    std::size_t const pdfs = pdf_uncertainty_pdfs_one.size() / 2;

    assert( scale_uncertainty_pdfs_one.size() == 2 * scales );
    assert( scale_uncertainty_pdfs_two.size() == 2 * scales );
    assert( matrix_elements.size() == 2 * scales );
    assert( alphas_factors.size() == 2 * scales );
    assert( pdf_uncertainty_pdfs_one.size() == 2 * pdfs );
    assert( pdf_uncertainty_pdfs_two.size() == 2 * pdfs );

    scale_results.clear();

    for (std::size_t i = scale_index; i != scale_index + scales; ++i)
    {
        scale_results.push_back(global_factor * alphas_factors.at(i) * convolute(
            scale_uncertainty_pdfs_one.at(i),
            scale_uncertainty_pdfs_two.at(i),
            matrix_elements.at(i),
            set
        ));
    }

    pdf_results.clear();

    for (std::size_t i = pdf_index; i != pdf_index + pdfs; ++i)
    {
        pdf_results.push_back(global_factor * alphas_factors.front() * convolute(
            pdf_uncertainty_pdfs_one.at(i),
            pdf_uncertainty_pdfs_two.at(i),
            matrix_elements.front(),
            set
        ));
    }
}

template <typename T, typename I>
inline void convolute_mes_with_pdfs(
    std::vector<neg_pos_results<T>>& scale_results,
    std::vector<neg_pos_results<T>>& pdf_results,
    std::vector<parton_array<T>> const& scale_uncertainty_pdfs_one,
    std::vector<parton_array<T>> const& scale_uncertainty_pdfs_two,
    std::vector<parton_array<T>> const& pdf_uncertainty_pdfs_one,
    std::vector<parton_array<T>> const& pdf_uncertainty_pdfs_two,
    std::vector<initial_state_map<T>> const& matrix_elements,
    initial_state_set set,
    std::vector<T> const& alphas_factors,
    T global_factor,
    cut_result_with_info<I> const& cut
) {
    std::size_t const scales = scale_uncertainty_pdfs_one.size();
    std::size_t const pdfs = pdf_uncertainty_pdfs_one.size();

    assert( scale_uncertainty_pdfs_one.size() == scales );
    assert( scale_uncertainty_pdfs_two.size() == scales );
    assert( matrix_elements.size() == scales );
    assert( alphas_factors.size() == scales );
    assert( pdf_uncertainty_pdfs_one.size() == pdfs );
    assert( pdf_uncertainty_pdfs_two.size() == pdfs );

    scale_results.clear();

    for (std::size_t i = 0; i != scales; ++i)
    {
        scale_results.push_back(convolute(
            scale_uncertainty_pdfs_one.at(i),
            scale_uncertainty_pdfs_two.at(i),
            matrix_elements.at(i),
            set,
            alphas_factors.at(i) * global_factor,
            cut
        ));
    }

    pdf_results.clear();

    for (std::size_t i = 0; i != pdfs; ++i)
    {
        pdf_results.push_back(convolute(
            pdf_uncertainty_pdfs_one.at(i),
            pdf_uncertainty_pdfs_two.at(i),
            matrix_elements.front(),
            set,
            global_factor,
            cut
        ));
    }
}

}

#endif
