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
        pdf_results_.reserve((pdfs_.count() == 1) ? 0 : pdfs_.count());

        if (!dynamic_scales_)
        {
            std::vector<T> no_point;
            std::vector<recombined_state> no_states;
            psp<T> no_psp{no_point, no_states, T(), psp_type::pos_rap};

            // a static scale must not depend on any phase or state information
            set_scales(no_psp, no_psp);

            results_.reserve(pos_.scales_.size());
            pos_.borns.resize(pos_.scales_.size());
            neg_.borns.resize(neg_.scales_.size());
        }

        pdfs_.register_partons(partons_in_initial_state_set(set));
    }

    T eval(
        std::vector<T> const& phase_space,
        luminosity_info<T> const& info,
        hep::projector<T>& projector
    ) override {
        // assumes that the reombination is independent of the frame
        recombiner_.recombine(phase_space, final_states_, pos_.ps, pos_.states);
        neg_.ps.assign(pos_.ps.begin(), pos_.ps.end());
        neg_.states.assign(pos_.states.begin(), pos_.states.end());

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
            results_.reserve(pos_.scales_.size());
            pos_.borns.resize(pos_.scales_.size());
            neg_.borns.resize(neg_.scales_.size());
        }

        T const factor = T(0.5) * hbarc2_ / info.energy_squared();

        pdfs_.eval(info.x1(), neg_.scales_, neg_.scale_pdf_x1, neg_.pdf_pdf_x1);
        pdfs_.eval(info.x2(), neg_.scales_, neg_.scale_pdf_x2, neg_.pdf_pdf_x2);
        pdfs_.eval(info.x1(), pos_.scales_, pos_.scale_pdf_x1, pos_.pdf_pdf_x1);
        pdfs_.eval(info.x2(), pos_.scales_, pos_.scale_pdf_x2, pos_.pdf_pdf_x2);

        assert( neg_.scale_pdf_x1.size() == neg_.scales_.size() );
        assert( neg_.scale_pdf_x2.size() == neg_.scales_.size() );
        assert( (pdfs_.count() == 1) || (neg_.pdf_pdf_x1.size() == pdfs_.count()) );
        assert( (pdfs_.count() == 1) || (neg_.pdf_pdf_x2.size() == pdfs_.count()) );
        assert( pos_.scale_pdf_x1.size() == pos_.scales_.size() );
        assert( pos_.scale_pdf_x2.size() == pos_.scales_.size() );
        assert( (pdfs_.count() == 1) || (pos_.pdf_pdf_x1.size() == pdfs_.count()) );
        assert( (pdfs_.count() == 1) || (pos_.pdf_pdf_x2.size() == pdfs_.count()) );

        for (auto& born : pos_.borns)
        {
            born.clear();
        }

        bool neg_pos_scales_different = true;
        std::size_t size;

        if (neg_pos_scales_different)
        {
            size = pos_.scales_.size();
            pos_.borns.resize(2 * size);
            pos_.scales_.insert(pos_.scales_.end(), neg_.scales_.begin(), neg_.scales_.end());
        }

        matrix_elements_.borns(phase_space, set_, pos_.scales_, pos_.borns);

        assert( pos_.borns.size() == pos_.scales_.size() );

        if (neg_pos_scales_different)
        {
            pos_.scales_.erase(std::prev(pos_.scales_.end(), size), pos_.scales_.end());
            neg_.borns.assign(std::prev(pos_.borns.end(), size), pos_.borns.end());
            pos_.borns.erase(std::prev(pos_.borns.end(), size), pos_.borns.end());
        }
        else
        {
            neg_.borns = pos_.borns;
        }

        T result = T();

        if (!pos_cutted)
        {
            convolute_mes_with_pdfs(
                results_,
                pdf_results_,
                pos_.scale_pdf_x1,
                pos_.scale_pdf_x2,
                pos_.pdf_pdf_x1,
                pos_.pdf_pdf_x2,
                pos_.borns,
                set_,
                pos_.factors,
                factor
            );

            distributions_(pos_psp, results_, pdf_results_, projector);

            result += results_.front();
        }

        if (!neg_cutted)
        {
            convolute_mes_with_pdfs(
                results_,
                pdf_results_,
                neg_.scale_pdf_x2,
                neg_.scale_pdf_x1,
                neg_.pdf_pdf_x2,
                neg_.pdf_pdf_x1,
                neg_.borns,
                set_,
                neg_.factors,
                factor
            );

            result += results_.front();

            distributions_(pos_psp, results_, pdf_results_, projector);
        }

        return result;
    }

protected:
    void set_scales(psp<T> const& neg_psp, psp<T> const& pos_psp)
    {
        using std::pow;

        neg_.scales_.clear();
        neg_.factors.clear();
        pos_.scales_.clear();
        pos_.factors.clear();
        scale_setter_(neg_psp, neg_.scales_);
        scale_setter_(pos_psp, pos_.scales_);
        pdfs_.eval_alphas(neg_.scales_, neg_.factors);
        pdfs_.eval_alphas(pos_.scales_, pos_.factors);

        T const central_alphas = pos_.factors.front();
        matrix_elements_.alphas(central_alphas);

        for (T& factor : neg_.factors)
        {
            T const alphas = factor;
            factor = pow(alphas / central_alphas, alphas_power_);
        }
        for (T& factor : pos_.factors)
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
        std::vector<parton_array<T>> scale_pdf_x1;
        std::vector<parton_array<T>> scale_pdf_x2;
        std::vector<parton_array<T>> pdf_pdf_x1;
        std::vector<parton_array<T>> pdf_pdf_x2;
        std::vector<scales<T>> scales_;
        std::vector<T> factors;
        std::vector<initial_state_map<T>> borns;
    } neg_, pos_;

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
