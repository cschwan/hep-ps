#ifndef HEP_PS_BORN_INTEGRAND_HPP
#define HEP_PS_BORN_INTEGRAND_HPP

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

#include "hep/mc/projector.hpp"

#include "hep/ps/convolute.hpp"
#include "hep/ps/final_state.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/ps_integrand.hpp"
#include "psp.hpp"
#include "psp_type.hpp"
#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scales.hpp"

#include "nonstd/span.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace hep
{

template <class T, class M, class C, class R, class P, class S, class D>
class born_integrand : public ps_integrand<T>
{
public:
    template <
        typename MatrixElements,
        typename Cuts,
        typename Recombiner,
        typename Pdfs,
        typename ScaleSetter,
        typename Distributions>
    born_integrand(
        MatrixElements&& matrix_elements,
        Cuts&& cuts,
        Recombiner&& recombiner,
        Pdfs&& pdfs,
        ScaleSetter&& scale_setter,
        Distributions&& distributions,
        initial_state_set set,
        T hbarc2
    )
        : matrix_elements_(std::forward<MatrixElements>(matrix_elements))
        , cuts_(std::forward<Cuts>(cuts))
        , recombiner_(std::forward<Recombiner>(recombiner))
        , pdfs_(std::forward<Pdfs>(pdfs))
        , scale_setter_(std::forward<ScaleSetter>(scale_setter))
        , distributions_(std::forward<Distributions>(distributions))
        , set_(set)
        , hbarc2_(hbarc2)
        , final_states_(matrix_elements_.final_states())
        , alphas_power_(matrix_elements_.alphas_power())
        , dynamic_scales_{scale_setter_.dynamic()}
    {
        std::size_t const fs = final_states_.size();

        neg_.ps.reserve(4 * (fs + 2));
        pos_.ps.reserve(4 * (fs + 2));
        neg_.states.reserve(fs);
        pos_.states.reserve(fs);

        std::size_t const scale_count = scale_setter_.count();
        std::size_t const pdf_count = pdfs_.count();

        results_.reserve(scale_count);
        pdf_results_.reserve((pdf_count == 1) ? 0 : pdf_count);
        borns_.resize(2 * scale_count);
        scales_.resize(2 * scale_count);

        if (!dynamic_scales_)
        {
            std::vector<T> no_point;
            std::vector<recombined_state> no_states;
            psp<T> no_psp{no_point, no_states, T(), psp_type::pos_rap};

            // static scales must not depend on any phase space or state information
            set_scales(no_psp, no_psp);
        }

        pdfs_.register_partons(partons_in_initial_state_set(set));
    }

    T eval(
        std::vector<T> const& phase_space,
        luminosity_info<T> const& info,
        hep::projector<T>& projector
    ) override {
        recombiner_.recombine(phase_space, final_states_, pos_.ps, pos_.states);

        // assumes that the reombination is independent of the frame
        neg_.ps = pos_.ps;
        neg_.states = pos_.states;

        psp<T> neg_psp{neg_.ps, neg_.states, info.rapidity_shift(), psp_type::neg_rap};
        psp<T> pos_psp{pos_.ps, pos_.states, info.rapidity_shift(), psp_type::pos_rap};

        bool const neg_cutted = cuts_.cut(neg_psp);
        bool const pos_cutted = cuts_.cut(pos_psp);

        if (neg_cutted && pos_cutted)
        {
            return T();
        }

        if (dynamic_scales_)
        {
            set_scales(neg_psp, pos_psp);
        }

        std::size_t const scales = scale_setter_.count();
        std::size_t const pdfs = (pdfs_.count() == 1) ? 0 : pdfs_.count();

        T const factor = T(0.5) * hbarc2_ / info.energy_squared();

        pdfs_.eval(info.x1(), scales, scales_, scale_pdf_x1_, pdf_pdf_x1_);
        pdfs_.eval(info.x2(), scales, scales_, scale_pdf_x2_, pdf_pdf_x2_);

        // for each of the two phase space point the must be a set of scales and PDFs
        assert( scales_.size() == 2 * scales );
        assert( scale_pdf_x1_.size() == 2 * scales );
        assert( scale_pdf_x2_.size() == 2 * scales );

        // for each CENTRAL scale there must be a set of uncertainty PDFs
        assert( pdf_pdf_x1_.size() == 2 * pdfs );
        assert( pdf_pdf_x2_.size() == 2 * pdfs );

        for (auto& born : borns_)
        {
            born.clear();
        }

        matrix_elements_.borns(phase_space, set_, nonstd::span<hep::scales<T>>{scales_}, borns_);

        // for each scale there must be a matrix element
        assert( borns_.size() == 2 * scales );

        T result = T();

        using span0 = nonstd::span<parton_array<T> const>;
        using span1 = nonstd::span<initial_state_map<T> const>;
        using span2 = nonstd::span<T const>;

        if (!neg_cutted)
        {
            convolute_mes_with_pdfs(
                results_,
                pdf_results_,
                span0{scale_pdf_x2_}.first(scales),
                span0{scale_pdf_x1_}.first(scales),
                span0{pdf_pdf_x2_}.first(pdfs),
                span0{pdf_pdf_x1_}.first(pdfs),
                span1{borns_}.first(scales),
                set_,
                span2{factors_}.first(scales),
                factor
            );

            result += results_.front();

            distributions_(neg_psp, results_, pdf_results_, projector);
        }

        if (!pos_cutted)
        {
            convolute_mes_with_pdfs(
                results_,
                pdf_results_,
                span0{scale_pdf_x1_}.last(scales),
                span0{scale_pdf_x2_}.last(scales),
                span0{pdf_pdf_x1_}.last(pdfs),
                span0{pdf_pdf_x2_}.last(pdfs),
                span1{borns_}.last(scales),
                set_,
                span2{factors_}.last(scales),
                factor
            );

            distributions_(pos_psp, results_, pdf_results_, projector);

            result += results_.front();
        }

        return result;
    }

protected:
    void set_scales(psp<T> const& neg_psp, psp<T> const& pos_psp)
    {
        using std::pow;
        using span = nonstd::span<scales<T>>;

        std::size_t const scales = scale_setter_.count();

        factors_.clear();

        scale_setter_.eval(neg_psp, span{scales_}.first(scales));
        scale_setter_.eval(pos_psp, span{scales_}.last(scales));

        pdfs_.eval_alphas(scales_, factors_);

        T const central_alphas = factors_.front();
        matrix_elements_.alphas(central_alphas);

        for (T& factor : factors_)
        {
            T const alphas = factor;
            factor = pow(alphas / central_alphas, alphas_power_);
        }
    }

private:
    M matrix_elements_;
    C cuts_;
    R recombiner_;
    P pdfs_;
    S scale_setter_;
    D distributions_;
    initial_state_set set_;
    T hbarc2_;

    struct
    {
        std::vector<T> ps;
        std::vector<recombined_state> states;
    } neg_, pos_;

    std::vector<parton_array<T>> scale_pdf_x1_;
    std::vector<parton_array<T>> scale_pdf_x2_;
    std::vector<parton_array<T>> pdf_pdf_x1_;
    std::vector<parton_array<T>> pdf_pdf_x2_;
    std::vector<scales<T>> scales_;
    std::vector<T> factors_;
    std::vector<initial_state_map<T>> borns_;
    std::vector<final_state> final_states_;
    std::vector<T> pdf_results_;
    std::vector<T> results_;
    T alphas_power_;
    bool dynamic_scales_;
};

template <class T, class M, class C, class R, class P, class S, class D>
using born_integrand_t = born_integrand<T, std::decay_t<M>, std::decay_t<C>, std::decay_t<R>,
    std::decay_t<P>, std::decay_t<S>, std::decay_t<D>>;

template <class T, class M, class C, class R, class P, class S, class D>
inline std::unique_ptr<ps_integrand<T>> make_born_integrand(
    M&& matrix_elements,
    C&& cuts,
    R&& recombiner,
    P&& pdfs,
    S&& scale_setter,
    D&& distributions,
    initial_state_set set,
    T hbarc2
) {
    return std::make_unique<born_integrand_t<T, M, C, R, P, S, D>>(
        std::forward<M>(matrix_elements),
        std::forward<C>(cuts),
        std::forward<R>(recombiner),
        std::forward<P>(pdfs),
        std::forward<S>(scale_setter),
        std::forward<D>(distributions),
        set,
        hbarc2
    );
}

}

#endif
