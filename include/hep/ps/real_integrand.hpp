#ifndef HEP_PS_REAL_INTEGRAND_HPP
#define HEP_PS_REAL_INTEGRAND_HPP

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
#include "hep/ps/non_zero_dipole.hpp"
#include "hep/ps/particle_type.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/ps_integrand.hpp"
#include "hep/ps/psp.hpp"
#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scales.hpp"

#include "nonstd/span.hpp"

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <memory>
#include <type_traits>
#include <vector>

namespace hep
{

template <class T, class M, class S, class C, class R, class P, class U, class D>
class real_integrand : public ps_integrand<T>
{
public:
    template <
        typename MatrixElements,
        typename Subtraction,
        typename Cuts,
        typename Recombiner,
        typename Pdfs,
        typename ScaleSetter,
        typename Distributions>
    real_integrand(
        MatrixElements&& matrix_elements,
        Subtraction&& subtraction,
        Cuts&& cuts,
        Recombiner&& recombiner,
        Pdfs&& pdfs,
        ScaleSetter&& scale_setter,
        Distributions&& distributions,
        initial_state_set set,
        T hbarc2,
        T alpha_min
    )
        : matrix_elements_(std::forward<MatrixElements>(matrix_elements))
        , subtraction_(std::forward<Subtraction>(subtraction))
        , cuts_(std::forward<Cuts>(cuts))
        , recombiner_(std::forward<Recombiner>(recombiner))
        , pdfs_(std::forward<Pdfs>(pdfs))
        , scale_setter_(std::forward<ScaleSetter>(scale_setter))
        , distributions_(std::forward<Distributions>(distributions))
        , set_(set)
        , hbarc2_(hbarc2)
        , alpha_min_(alpha_min)
        , phase_space_indices_()
        , phase_space_sizes_{1}
        , final_states_real_(matrix_elements_.final_states_real())
        , final_states_dipole_(matrix_elements_.final_states())
        , dipoles_(matrix_elements_.dipoles())
        , alphas_power_(matrix_elements_.alphas_power())
    {
        std::size_t const fs = final_states_real_.size();

        recombined_ps_.reserve(4 * (fs + 2));
        recombined_states_.reserve(fs);
        recombined_dipole_states_.reserve(fs - 1);
        pdf_results_.reserve((pdfs_.count() == 1) ? 0 : pdfs_.count());

        non_zero_dipoles_.reserve(dipoles_.size());

        pdfs_.register_partons(partons_in_initial_state_set(set));

        for (std::size_t i = 0; i != dipoles_.size() - 1; ++i)
        {
            auto const& a = dipoles_.at(i);

            for (std::size_t j = i + 1; j != dipoles_.size(); ++j)
            {
                auto const& b = dipoles_.at(j);

                if (subtraction_.same_mapping(a, b))
                {
                    std::rotate(
                        std::next(dipoles_.begin(), i + 1),
                        std::next(dipoles_.begin(), j),
                        std::next(dipoles_.begin(), j + 1)
                    );

                    break;
                }
            }
        }

        for (std::size_t i = 1; i != dipoles_.size(); ++i)
        {
            auto const& a = dipoles_.at(i - 1);
            auto const& b = dipoles_.at(i);

            if (subtraction_.same_mapping(a, b))
            {
                ++phase_space_sizes_.back();
            }
            else
            {
                phase_space_sizes_.push_back(1);
            }
        }

        // check if `cuts` can hold the cut results for dipole phase space
        assert( neg_.pass_cut.size() >= phase_space_sizes_.size() );
        assert( neg_.pass_cut.size() >= phase_space_sizes_.size() );

        std::size_t const scales = scale_setter_.count();
        scales_.resize(2 * scales * (phase_space_sizes_.size() + 1));
        factors_.reserve(scales_.size());

        if (!scale_setter_.dynamic())
        {
            std::vector<T> no_point;
            std::vector<recombined_state> no_states;
            psp<T> no_psp{no_point, no_states, T(), psp_type::pos_rap};

            // static scales must not depend on any phase space or state information
            using span0 = nonstd::span<hep::scales<T>>;

            scale_setter_.eval(no_psp, span0{scales_}.first(scales));

            for (std::size_t i = 1; i != scales_.size() / scales; ++i)
            {
                std::copy(scales_.begin(), std::next(scales_.begin(), scales),
                    std::next(scales_.begin(), i * scales));
            }

            pdfs_.eval_alphas(scales_, factors_);

            T const central_alphas = factors_.front();
            matrix_elements_.alphas(central_alphas);

            using std::pow;

            for (T& factor : factors_)
            {
                T const alphas = factor;
                factor = pow(alphas / central_alphas, alphas_power_);
            }
        }

        results_.reserve(scales);
        me_.resize(2 * scales);
        me_tmp_.resize(2 * scales);

        phase_space_indices_.reserve(phase_space_sizes_.size());
        dipole_phase_spaces_.resize(phase_space_sizes_.size());
        neg_.ps.resize(phase_space_sizes_.size());
        pos_.ps.resize(phase_space_sizes_.size());
        neg_.states.resize(phase_space_sizes_.size());
        pos_.states.resize(phase_space_sizes_.size());
        neg_.real_ps.resize(4 * (fs + 2));
        neg_.real_states.resize(fs);
        pos_.real_ps.resize(4 * (fs + 2));
        pos_.real_states.resize(fs);

        for (auto& ps : dipole_phase_spaces_)
        {
            ps.resize(4 * (fs + 1));
        }

        for (auto& ps : neg_.ps)
        {
            ps.resize(4 * (fs + 1));
        }

        for (auto& ps : pos_.ps)
        {
            ps.resize(4 * (fs + 1));
        }

        for (auto& states : neg_.states)
        {
            states.reserve(fs - 1);
        }

        for (auto& states : pos_.states)
        {
            states.reserve(fs - 1);
        }
    }

    T eval(
        std::vector<T> const& real_phase_space,
        luminosity_info<T> const& info,
        hep::projector<T>& projector
    ) override {
        non_zero_dipoles_.clear();
        phase_space_indices_.clear();

        neg_.pass_cut.reset();
        pos_.pass_cut.reset();

        std::size_t dipole_index = 0;
        std::size_t phase_space_index = 0;

        // STEP 1: generate the different phase spaces for all dipoles and check whether the dipoles
        // are cut away, the whole point is tech-cutted, or whether we have to evaluate it

        for (std::size_t phase_space_size : phase_space_sizes_)
        {
            auto& cms_ps = dipole_phase_spaces_.at(phase_space_index);
            auto& neg_ps = neg_.ps.at(phase_space_index);
            auto& pos_ps = pos_.ps.at(phase_space_index);
            auto& neg_states = neg_.states.at(phase_space_index);
            auto& pos_states = pos_.states.at(phase_space_index);

            bool neg_cutted = true;
            bool pos_cutted = true;

            for (std::size_t i = 0; i != phase_space_size; ++i)
            {
                auto const& dipole = dipoles_.at(dipole_index++);

                // map the real phase space on the dipole phase space
                auto const inv = subtraction_.map_phase_space(real_phase_space, cms_ps, dipole);

                if (inv.alpha < alpha_min_)
                {
                    // if we apply a technical cut, completely throw away the point
                    return T();
                }

                if (i == 0)
                {
                    recombiner_.recombine(cms_ps, final_states_dipole_, pos_ps, pos_states);

                    // assumes that the recombination is independent of the frame
                    neg_ps = pos_ps;
                    neg_states = pos_states;

                    psp<T> neg_psp{neg_ps, neg_states, info.rapidity_shift(), psp_type::neg_rap};
                    psp<T> pos_psp{pos_ps, pos_states, info.rapidity_shift(), psp_type::pos_rap};

                    neg_cutted = cuts_.cut(neg_psp);
                    pos_cutted = cuts_.cut(pos_psp);

                    ++phase_space_index;

                    if (neg_cutted && pos_cutted)
                    {
                        // if the dipole is cutted, we still have to check the technical cuts
                        continue;
                    }

                    neg_.pass_cut.set(phase_space_index - 1, !neg_cutted);
                    pos_.pass_cut.set(phase_space_index - 1, !pos_cutted);

                    phase_space_indices_.push_back(phase_space_index - 1);
                }

                if (!neg_cutted || !pos_cutted)
                {
                    non_zero_dipoles_.emplace_back(inv, dipole);
                }
            }
        }

        // STEP 2: Check if the real matrix element has to evaluated and check if this event is
        // nonzero

        recombiner_.recombine(real_phase_space, final_states_real_, pos_.real_ps, pos_.real_states);

        // assumes that the recombination is independent of the frame
        neg_.real_ps = pos_.real_ps;
        neg_.real_states = pos_.real_states;

        psp<T> neg_psp{neg_.real_ps, neg_.real_states, info.rapidity_shift(), psp_type::neg_rap};
        psp<T> pos_psp{pos_.real_ps, pos_.real_states, info.rapidity_shift(), psp_type::pos_rap};

        bool const neg_cutted = cuts_.cut(neg_psp);
        bool const pos_cutted = cuts_.cut(pos_psp);

        if (phase_space_indices_.empty() && neg_cutted && pos_cutted)
        {
            return T();
        }

        // STEP 3: Calculate the scales, alphas, and PDFs for all dipoles

        std::size_t const scales = scale_setter_.count();
        std::size_t const pdfs = (pdfs_.count() == 1) ? 0 : pdfs_.count();

        using span4 = nonstd::span<hep::scales<T> const>;
        auto scale_span = span4{scales_}.first(2 * scales * (phase_space_indices_.size() + 1));

        if (scale_setter_.dynamic())
        {
            factors_.clear();

            using span0 = nonstd::span<hep::scales<T>>;

            // evaluate scales for the dipoles, but only if they passed the cuts
            for (std::size_t i = 0; i != phase_space_indices_.size(); ++i)
            {
                auto const phase_space_index = phase_space_indices_.at(i);
                auto neg_span = span0{scales_}.subspan(scales * (2 * i + 0), scales);
                auto pos_span = span0{scales_}.subspan(scales * (2 * i + 1), scales);

                if (neg_.pass_cut.test(phase_space_index))
                {
                    auto& neg_ps = neg_.ps.at(phase_space_index);
                    auto& neg_states = neg_.states.at(phase_space_index);
                    psp<T> const neg_psp{neg_ps, neg_states, info.rapidity_shift(),
                        psp_type::neg_rap};

                    scale_setter_.eval(neg_psp, neg_span);
                }

                if (pos_.pass_cut.test(phase_space_index))
                {
                    auto& pos_ps = pos_.ps.at(phase_space_index);
                    auto& pos_states = pos_.states.at(phase_space_index);
                    psp<T> const pos_psp{pos_ps, pos_states, info.rapidity_shift(),
                        psp_type::pos_rap};

                    scale_setter_.eval(pos_psp, pos_span);

                    if (!neg_.pass_cut.test(phase_space_index))
                    {
                        std::copy(pos_span.begin(), pos_span.end(), neg_span.begin());
                    }
                }
                else
                {
                    std::copy(neg_span.begin(), neg_span.end(), pos_span.begin());
                }
            }

            auto neg_span = span0{scales_}.subspan(scales * (2 * phase_space_indices_.size() + 0),
                scales);
            auto pos_span = span0{scales_}.subspan(scales * (2 * phase_space_indices_.size() + 1),
                scales);

            // evaluate scales for the real matrix element

            if (!neg_cutted || !pos_cutted)
            {
                if (!neg_cutted)
                {
                    auto& neg_ps = neg_.real_ps;
                    auto& neg_states = neg_.real_states;
                    psp<T> const neg_psp{neg_ps, neg_states, info.rapidity_shift(),
                        psp_type::neg_rap};

                    scale_setter_.eval(neg_psp, neg_span);
                }

                if (!pos_cutted)
                {
                    auto& pos_ps = pos_.real_ps;
                    auto& pos_states = pos_.real_states;
                    psp<T> const pos_psp{pos_ps, pos_states, info.rapidity_shift(),
                        psp_type::pos_rap};

                    scale_setter_.eval(pos_psp, pos_span);

                    if (neg_cutted)
                    {
                        std::copy(pos_span.begin(), pos_span.end(), neg_span.begin());
                    }
                }
                else
                {
                    std::copy(neg_span.begin(), neg_span.end(), pos_span.begin());
                }
            }
            else
            {
                auto span = span0{scales_}.subspan(scales * 2 * (phase_space_indices_.size() - 1),
                    scales);

                std::copy(span.begin(), span.end(), neg_span.begin());
                std::copy(span.begin(), span.end(), pos_span.begin());
            }

            pdfs_.eval_alphas(scale_span, factors_);

            T const central_alphas = factors_.front();
            matrix_elements_.alphas(central_alphas);

            using std::pow;

            for (T& factor : factors_)
            {
                T const alphas = factor;
                factor = pow(alphas / central_alphas, alphas_power_);
            }
        }

        pdfs_.eval(info.x1(), scales, scale_span, pdfsx1_, pdf_pdfsx1_);
        pdfs_.eval(info.x2(), scales, scale_span, pdfsx2_, pdf_pdfsx2_);

        assert( pdfsx1_.size() == 2 * (phase_space_indices_.size() + 1) * scales );
        assert( pdfsx2_.size() == 2 * (phase_space_indices_.size() + 1) * scales );
        assert( pdf_pdfsx1_.size() == 2 * (phase_space_indices_.size() + 1) * pdfs );
        assert( pdf_pdfsx2_.size() == 2 * (phase_space_indices_.size() + 1) * pdfs );

        T const factor = T(0.5) * hbarc2_ / info.energy_squared();

        T dipoles = T();

        std::size_t non_zero_dipole_index = 0;

        using span1 = nonstd::span<parton_array<T> const>;
        using span2 = nonstd::span<initial_state_map<T> const>;
        using span3 = nonstd::span<T const>;

        for (std::size_t i = 0; i != phase_space_indices_.size(); ++i)
        {
            auto const phase_space_index = phase_space_indices_.at(i);
            auto const& cms_ps = dipole_phase_spaces_.at(phase_space_index);
            auto& neg_ps = neg_.ps.at(phase_space_index);
            auto& pos_ps = pos_.ps.at(phase_space_index);
            auto& neg_states = neg_.states.at(phase_space_index);
            auto& pos_states = pos_.states.at(phase_space_index);

            psp<T> const neg_psp{neg_ps, neg_states, info.rapidity_shift(), psp_type::neg_rap};
            psp<T> const pos_psp{pos_ps, pos_states, info.rapidity_shift(), psp_type::pos_rap};

            for (std::size_t j = 0; j != phase_space_sizes_.at(phase_space_index); ++j)
            {
                auto const& non_zero_dipole = non_zero_dipoles_.at(non_zero_dipole_index++);
                auto const& invariants = non_zero_dipole.invariants();
                auto const& dipole = non_zero_dipole.dipole();

                for (auto& me : me_)
                {
                    me.clear();
                }

                std::size_t const scale_index = scales * 2 * i;
                std::size_t const pdf_index = pdfs * 2 * i;

                T function = T(1.0);

                if ((dipole.emitter_type() == particle_type::fermion) !=
                    (dipole.unresolved_type() == particle_type::fermion))
                {
                    matrix_elements_.dipole_me(dipole, cms_ps, set_,
                        span4{scales_}.subspan(scale_index, 2 * scales), me_);
                    function = -subtraction_.fermion_function(dipole, invariants);
                }
                else
                {
                    auto const correlator = subtraction_.boson_function(dipole, invariants,
                        real_phase_space);

                    for (auto& me : me_tmp_)
                    {
                        me.clear();
                    }

                    matrix_elements_.dipole_sc(
                        dipole,
                        cms_ps,
                        correlator.p,
                        set_,
                        span4{scales_}.subspan(scale_index, 2 * scales),
                        me_,
                        me_tmp_
                    );

                    for (std::size_t k = 0; k != me_.size(); ++k)
                    {
                        for (std::size_t l = 0; l != me_[k].size(); ++l)
                        {
                            assert( me_[k][l].first == me_tmp_[k][l].first );

                            me_[k][l].second *= -correlator.a;
                            me_[k][l].second -= me_tmp_[k][l].second * correlator.b;
                        }
                    }
                }

                if (neg_.pass_cut.test(phase_space_index))
                {
                    convolute_mes_with_pdfs(
                        results_,
                        pdf_results_,
                        span1{pdfsx2_}.subspan(2 * scales * i, scales),
                        span1{pdfsx1_}.subspan(2 * scales * i, scales),
                        span1{pdf_pdfsx2_}.subspan(pdf_index, pdfs),
                        span1{pdf_pdfsx1_}.subspan(pdf_index, pdfs),
                        span2{me_}.first(scales),
                        set_,
                        span3{factors_}.subspan(2 * scales * i, scales),
                        function * factor
                    );

                    distributions_(neg_psp, results_, pdf_results_, projector);

                    dipoles += results_.front();
                }

                if (pos_.pass_cut.test(phase_space_index))
                {
                    convolute_mes_with_pdfs(
                        results_,
                        pdf_results_,
                        span1{pdfsx1_}.subspan(2 * scales * i + scales, scales),
                        span1{pdfsx2_}.subspan(2 * scales * i + scales, scales),
                        span1{pdf_pdfsx1_}.subspan(pdf_index + pdfs, pdfs),
                        span1{pdf_pdfsx2_}.subspan(pdf_index + pdfs, pdfs),
                        span2{me_}.last(scales),
                        set_,
                        span3{factors_}.subspan(2 * scales * i + scales, scales),
                        function * factor
                    );

                    distributions_(pos_psp, results_, pdf_results_, projector);

                    dipoles += results_.front();
                }
            }
        }

        T real = T();

        if (!neg_cutted || !pos_cutted)
        {
            for (auto& me : me_)
            {
                me.clear();
            }

            std::size_t const scale_index = scales * 2 * phase_space_indices_.size();
            std::size_t const pdf_index = pdfs * 2 * phase_space_indices_.size();

            matrix_elements_.reals(
                real_phase_space,
                set_,
                span4{scales_}.subspan(scale_index, 2 * scales),
                me_
            );

            if (!neg_cutted)
            {
                convolute_mes_with_pdfs(
                    results_,
                    pdf_results_,
                    span1{pdfsx2_}.subspan(scale_index, scales),
                    span1{pdfsx1_}.subspan(scale_index, scales),
                    span1{pdf_pdfsx2_}.subspan(pdf_index, pdfs),
                    span1{pdf_pdfsx1_}.subspan(pdf_index, pdfs),
                    span2{me_}.first(scales),
                    set_,
                    span3{factors_}.subspan(scale_index, scales),
                    factor
                );

                distributions_(neg_psp, results_, pdf_results_, projector);

                real += results_.front();
            }

            if (!pos_cutted)
            {
                convolute_mes_with_pdfs(
                    results_,
                    pdf_results_,
                    span1{pdfsx1_}.subspan(scale_index + scales, scales),
                    span1{pdfsx2_}.subspan(scale_index + scales, scales),
                    span1{pdf_pdfsx1_}.subspan(pdf_index + pdfs, pdfs),
                    span1{pdf_pdfsx2_}.subspan(pdf_index + pdfs, pdfs),
                    span2{me_}.last(scales),
                    set_,
                    span3{factors_}.subspan(scale_index + scales, scales),
                    factor
                );

                distributions_(pos_psp, results_, pdf_results_, projector);

                real += results_.front();
            }
        }

        return real + dipoles;
    }

private:
    M matrix_elements_;
    S subtraction_;
    C cuts_;
    R recombiner_;
    P pdfs_;
    U scale_setter_;
    D distributions_;
    initial_state_set set_;
    T hbarc2_;
    T alpha_min_;

    struct
    {
        std::vector<std::vector<T>> ps;
        std::vector<T> real_ps;
        std::vector<recombined_state> real_states;
        std::vector<std::vector<recombined_state>> states;
        std::vector<std::vector<scales<T>>> scales_;
        std::vector<scales<T>> real_scales;
        std::bitset<64> pass_cut;
    } neg_, pos_;

    std::vector<T> recombined_ps_;
    std::vector<std::vector<T>> dipole_phase_spaces_;
    std::vector<std::size_t> phase_space_indices_;
    std::vector<std::size_t> phase_space_sizes_;
    std::vector<recombined_state> recombined_states_;
    std::vector<recombined_state> recombined_dipole_states_;
    std::vector<recombined_state> dipole_recombined_states_;
    std::vector<final_state> final_states_real_;
    std::vector<final_state> final_states_dipole_;
    std::vector<parton_array<T>> pdfsx1_;
    std::vector<parton_array<T>> pdfsx2_;
    std::vector<parton_array<T>> pdf_pdfsx1_;
    std::vector<parton_array<T>> pdf_pdfsx2_;
    std::vector<T> pdf_results_;
    std::vector<T> results_;
    std::vector<scales<T>> scales_;
    std::vector<T> factors_;
    std::vector<initial_state_map<T>> me_;
    std::vector<initial_state_map<T>> me_tmp_;
    std::vector<non_zero_dipole<T>> non_zero_dipoles_;
    std::vector<dipole> dipoles_;
    T alphas_power_;
};

template <class T, class M, class S, class C, class R, class P, class U, class D>
using real_integrand_type = real_integrand<T, std::decay_t<M>, std::decay_t<S>, std::decay_t<C>,
    std::decay_t<R>, std::decay_t<P>, std::decay_t<U>, std::decay_t<D>>;

template <class T, class M, class S, class C, class R, class P, class U, class D>
inline std::unique_ptr<ps_integrand<T>> make_real_integrand(
    M&& matrix_elements,
    S&& subtraction,
    C&& cuts,
    R&& recombiner,
    P&& pdfs,
    U&& scale_setter,
    D&& distributions,
    initial_state_set set,
    T hbarc2,
    T alpha_min = T()
) {
    return std::make_unique<real_integrand_type<T, M, S, C, R, P, U, D>>(
        std::forward<M>(matrix_elements),
        std::forward<S>(subtraction),
        std::forward<C>(cuts),
        std::forward<R>(recombiner),
        std::forward<P>(pdfs),
        std::forward<U>(scale_setter),
        std::forward<D>(distributions),
        set,
        hbarc2,
        alpha_min
    );
}

}

#endif
