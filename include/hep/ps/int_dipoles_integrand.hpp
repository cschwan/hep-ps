#ifndef HEP_PS_INT_DIPOLES_INTEGRAND_HPP
#define HEP_PS_INT_DIPOLES_INTEGRAND_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2019  Christopher Schwan
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

#include "hep/ps/ab_terms.hpp"
#include "hep/ps/convolute.hpp"
#include "hep/ps/finite_parts.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/int_dipole.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/pdg_functions.hpp"
#include "hep/ps/ps_integrand.hpp"
#include "hep/ps/psp.hpp"
#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scales.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace hep
{

template <class T, class M, class S, class C, class R, class P, class U, class D>
class int_dipoles_integrand : public ps_integrand<T>
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
    int_dipoles_integrand(
        MatrixElements&& matrix_elements,
        Subtraction&& subtraction,
        Cuts&& cuts,
        Recombiner&& recombiner,
        Pdfs&& pdfs,
        ScaleSetter&& scale_setter,
        Distributions&& distributions,
        initial_state_set set,
        T hbarc2,
        finite_parts parts
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
        , parts_{parts}
        , final_states_(matrix_elements_.final_states())
        , insertion_terms_(matrix_elements_.dipoles())
        , alphas_power_(matrix_elements_.alphas_power())
    {
        std::size_t const fs = final_states_.size();

        neg_.ps.reserve(4 * (fs + 2));
        pos_.ps.reserve(4 * (fs + 2));
        neg_.states.reserve(fs);
        pos_.states.reserve(fs);

        std::size_t const scales = scale_setter_.count();

        neg_.results.reserve(scales);
        pos_.results.reserve(scales);
        neg_.pdf_results.reserve((pdfs_.count() == 1) ? 0 : pdfs_.count());
        pos_.pdf_results.reserve((pdfs_.count() == 1) ? 0 : pdfs_.count());
        scales_.resize(2 * scales);
        ab_terms_.reserve(scales);
        terms2_.reserve(scales);
        corr_me_.resize(insertion_terms_.size());

        if (!scale_setter_.dynamic())
        {
            std::vector<T> no_point;
            std::vector<recombined_state> no_states;
            psp<T> no_psp{no_point, no_states, T(), psp_type::pos_rap};

            // static scales must not depend on any phase space or state information
            set_scales(no_psp, false, no_psp, false);
        }

        pdfs_.register_partons(partons_in_initial_state_set(set));

        for (auto const state : set)
        {
            me_set_.add(state);

            auto const one = state_parton_one(state);
            auto const two = state_parton_two(state);

            switch (parton_type_of(one))
            {
            case parton_type::quark:
            case parton_type::anti_quark:
                me_set_.add(partons_to_initial_state(parton::photon, two));
                me_set_.add(partons_to_initial_state(parton::gluon, two));

                break;

            case parton_type::gluon_:
            case parton_type::photon_:
                for (auto p : parton_list())
                {
                    if ((p == parton::gluon) || (p == parton::photon))
                    {
                        continue;
                    }

                    me_set_.add(partons_to_initial_state(p, two));
                }

                break;

            default:
                assert( false );
            }

            switch (parton_type_of(two))
            {
            case parton_type::quark:
            case parton_type::anti_quark:
                me_set_.add(partons_to_initial_state(one, parton::photon));
                me_set_.add(partons_to_initial_state(one, parton::gluon));

                break;

            case parton_type::gluon_:
            case parton_type::photon_:
                for (auto p : parton_list())
                {
                    if ((p == parton::gluon) || (p == parton::photon))
                    {
                        continue;
                    }

                    me_set_.add(partons_to_initial_state(one, p));
                }

                break;

            default:
                assert( false );
            }
        }
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

        if (scale_setter_.dynamic())
        {
            set_scales(neg_psp, neg_cutted, pos_psp, pos_cutted);
        }

        std::size_t const scales = scale_setter_.count();
        std::size_t const pdfs = (pdfs_.count() == 1) ? 0 : pdfs_.count();

        neg_.results.clear();
        neg_.results.resize(scales);
        neg_.pdf_results.clear();
        neg_.pdf_results.resize(pdfs);
        pos_.results.clear();
        pos_.results.resize(scales);
        pos_.pdf_results.clear();
        pos_.pdf_results.resize(pdfs);

        T const x = phase_space.back();
        T const x1p = info.x1() / (info.x1() * (T(1.0) - x) + x);
        T const x2p = info.x2() / (info.x2() * (T(1.0) - x) + x);

        using span0 = nonstd::span<hep::scales<T> const>;

        auto const all_scales = span0{scales_}.subspan(
            scales * (neg_cutted ? 1 : 0),
            scales * ((!neg_cutted && !pos_cutted) ? 2 : 1)
        );

        pdfs_.eval(info.x1(), scales, all_scales, pdfsa_[0], pdf_pdfsa_[0]);
        pdfs_.eval(info.x2(), scales, all_scales, pdfsa_[1], pdf_pdfsa_[1]);
        pdfs_.eval(x1p, scales, all_scales, pdfsb_[0], pdf_pdfsb_[0]);
        pdfs_.eval(x2p, scales, all_scales, pdfsb_[1], pdf_pdfsb_[1]);

        for (auto& me : corr_me_)
        {
            me.clear();
        }

        // TODO: for the time being we assume that the matrix elements are independent of scales
        matrix_elements_.correlated_me(phase_space, me_set_, corr_me_);
        auto const factor = T(0.5) * hbarc2_ / info.energy_squared();

        T const eta[] = { info.x1(), info.x2() };
        T const xprime[] = {
            eta[0] * (T(1.0) - x) + x,
            eta[1] * (T(1.0) - x) + x
        };

        auto const neg_off_s = 0;
        auto const pos_off_s = neg_cutted ? 0 : scales;
        auto const neg_off_p = 0;
        auto const pos_off_p = neg_cutted ? 0 : pdfs;

        // loop over all Born, FI, IF, II, and FF
        for (std::size_t index = 0; index != insertion_terms_.size(); ++index)
        {
            auto const& me = corr_me_.at(index);
            auto const& term = insertion_terms_.at(index);

            // TODO: change interface of `correlated_me` based on the following assumption
            assert( me.size() == 1 );

            if ((parts_ != finite_parts::insertion_term2) &&
                (term.type() != insertion_term_type::final_final))
            {
                std::size_t const i = term.initial_particle();

                if (!neg_cutted)
                {
                    subtraction_.insertion_terms(term, span0{scales_}.first(scales), phase_space,
                        xprime[1 - i], eta[1 - i], ab_terms_);

                    assert( ab_terms_.size() == scales );

                    for (std::size_t j = 0; j != scales; ++j)
                    {
                        auto const eff_pdf = effective_pdf(
                            me.at(0).first,
                            ab_terms_.at(j),
                            term,
                            xprime[1 - i],
                            eta[1 - i],
                            pdfsa_[1 - i].at(j + neg_off_s),
                            pdfsb_[1 - i].at(j + neg_off_s)
                        );

                        neg_.results.at(j) += factors_.at(j) * factor * convolute(
                            (i == 0) ? eff_pdf                     : pdfsa_[1].at(j + neg_off_s),
                            (i == 0) ? pdfsa_[0].at(j + neg_off_s) : eff_pdf,
                            me,
                            me_set_
                        );
                    }

                    for (std::size_t j = 0; j != pdfs; ++j)
                    {
                        auto const eff_pdf = effective_pdf(
                            me.at(0).first,
                            ab_terms_.front(),
                            term,
                            xprime[1 - i],
                            eta[1 - i],
                            pdf_pdfsa_[1 - i].at(j + neg_off_p),
                            pdf_pdfsb_[1 - i].at(j + neg_off_p)
                        );

                        neg_.pdf_results.at(j) += factors_.at(0) * factor * convolute(
                            (i == 0) ? eff_pdf                         : pdf_pdfsa_[1].at(j + neg_off_p),
                            (i == 0) ? pdf_pdfsa_[0].at(j + neg_off_p) : eff_pdf,
                            me,
                            me_set_
                        );
                    }
                }

                if (!pos_cutted)
                {
                    subtraction_.insertion_terms(term, span0{scales_}.last(scales), phase_space,
                        xprime[i], eta[i], ab_terms_);

                    assert( ab_terms_.size() == scales );

                    for (std::size_t j = 0; j != scales; ++j)
                    {
                        auto const eff_pdf = effective_pdf(
                            me.at(0).first,
                            ab_terms_.at(j),
                            term,
                            xprime[i],
                            eta[i],
                            pdfsa_[i].at(j + pos_off_s),
                            pdfsb_[i].at(j + pos_off_s)
                        );

                        pos_.results.at(j) += factors_.at(j + neg_off_s) * factor * convolute(
                            (i == 0) ? eff_pdf                     : pdfsa_[0].at(j + pos_off_s),
                            (i == 0) ? pdfsa_[1].at(j + pos_off_s) : eff_pdf,
                            me,
                            me_set_
                        );
                    }

                    for (std::size_t j = 0; j != pdfs; ++j)
                    {
                        auto const eff_pdf = effective_pdf(
                            me.at(0).first,
                            ab_terms_.front(),
                            term,
                            xprime[i],
                            eta[i],
                            pdf_pdfsa_[i].at(j + pos_off_p),
                            pdf_pdfsb_[i].at(j + pos_off_p)
                        );

                        pos_.pdf_results.at(j) += factors_.at(0 + neg_off_s) * factor * convolute(
                            (i == 0) ? eff_pdf                         : pdf_pdfsa_[0].at(j + pos_off_p),
                            (i == 0) ? pdf_pdfsa_[1].at(j + pos_off_p) : eff_pdf,
                            me,
                            me_set_
                        );
                    }
                }
            }

            if ((parts_ != finite_parts::insertion_term) &&
                (term.type() != insertion_term_type::born))
            {
                subtraction_.insertion_terms2(term, all_scales, phase_space, terms2_);

                if (!neg_cutted)
                {
                    for (std::size_t i = 0; i != scales; ++i)
                    {
                        auto const f = factors_.at(i + neg_off_s) * factor * terms2_.at(i + neg_off_s);

                        neg_.results.at(i) += f * convolute(
                            pdfsa_[1].at(i + neg_off_s),
                            pdfsa_[0].at(i + neg_off_s),
                            me,
                            set_
                        );
                    }

                    for (std::size_t i = 0; i != pdfs; ++i)
                    {
                        auto const f = factors_.at(0 + neg_off_s) * factor * terms2_.at(0 + neg_off_s);

                        neg_.pdf_results.at(i) += f * convolute(
                            pdf_pdfsa_[1].at(i + neg_off_p),
                            pdf_pdfsa_[0].at(i + neg_off_p),
                            me,
                            set_
                        );
                    }
                }

                if (!pos_cutted)
                {
                    for (std::size_t i = 0; i != scales; ++i)
                    {
                        auto const f = factors_.at(i + pos_off_s) * factor * terms2_.at(i + pos_off_s);

                        pos_.results.at(i) += f * convolute(
                            pdfsa_[0].at(i + pos_off_s),
                            pdfsa_[1].at(i + pos_off_s),
                            me,
                            set_
                        );
                    }

                    for (std::size_t i = 0; i != pdfs; ++i)
                    {
                        auto const f = factors_.at(0 + pos_off_s) * factor * terms2_.at(0 + pos_off_s);

                        pos_.pdf_results.at(i) += f * convolute(
                            pdf_pdfsa_[0].at(i + pos_off_p),
                            pdf_pdfsa_[1].at(i + pos_off_p),
                            me,
                            set_
                        );
                    }
                }
            }
        }

        T result = T();

        if (!neg_cutted)
        {
            result += neg_.results.front();
            distributions_(neg_psp, neg_.results, neg_.pdf_results, projector);
        }

        if (!pos_cutted)
        {
            result += pos_.results.front();
            distributions_(pos_psp, pos_.results, pos_.pdf_results, projector);
        }

        return result;
    }

protected:
    parton_array<T> effective_pdf(
        initial_state me_in_state,
        ab_term<T> const& ab,
        int_dipole const& term,
        T xprime,
        T eta,
        parton_array<T> const& pdfa,
        parton_array<T> const& pdfb
    ) {
        parton_array<T> pdf;

        if (term.type() != insertion_term_type::final_initial)
        {
            auto a = pdg_id_to_parton(term.vertex().internal());
            auto ap = pdg_id_to_parton(term.vertex().external());

            T const convolute_factor =
                (state_parton_one(me_in_state) == state_parton_two(me_in_state)) ? T(0.5) : T(1.0);

            T actual_factor = T(1.0);

            if (((term.initial_particle() == 0) && (state_parton_one(me_in_state) == ap)) ||
                ((term.initial_particle() == 1) && (state_parton_two(me_in_state) == ap)))
            {
                actual_factor = T(0.5);
            }

            T const correction = actual_factor / convolute_factor;

            pdf[a]  = correction * pdfb[ap] * ab.a * (T(1.0) - eta) / xprime;
            pdf[a] += correction * pdfa[ap] * ab.b;
        }
        else
        {
            for (auto const a : parton_list())
            {
                pdf[a]  = pdfb[a] * ab.a * (T(1.0) - eta) / xprime;
                pdf[a] += pdfa[a] * ab.b;
            }
        }

        return pdf;
    }

    void set_scales(psp<T> const& neg_psp, bool neg_cutted, psp<T> const& pos_psp, bool pos_cutted)
    {
        using std::pow;
        using span = nonstd::span<scales<T>>;

        std::size_t const scales = scale_setter_.count();

        factors_.clear();

        if (!neg_cutted)
        {
            scale_setter_.eval(neg_psp, span{scales_}.first(scales));
        }

        if (!pos_cutted)
        {
            scale_setter_.eval(pos_psp, span{scales_}.last(scales));
        }

        auto const offset = neg_cutted ? scales : 0;
        auto const size = (!neg_cutted && !pos_cutted) ? 2 * scales : scales;

        pdfs_.eval_alphas(span{scales_}.subspan(offset, size), factors_);

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
    S subtraction_;
    C cuts_;
    R recombiner_;
    P pdfs_;
    U scale_setter_;
    D distributions_;
    initial_state_set set_;
    initial_state_set me_set_;
    T hbarc2_;
    finite_parts parts_;

    struct
    {
        std::vector<T> ps;
        std::vector<recombined_state> states;
        std::vector<T> results;
        std::vector<T> pdf_results;
    } neg_, pos_;

    std::vector<final_state> final_states_;
    std::array<std::vector<parton_array<T>>, 2> pdfsa_;
    std::array<std::vector<parton_array<T>>, 2> pdfsb_;
    std::array<std::vector<parton_array<T>>, 2> pdf_pdfsa_;
    std::array<std::vector<parton_array<T>>, 2> pdf_pdfsb_;
    std::vector<scales<T>> scales_;
    std::vector<T> factors_;
    std::vector<initial_state_map<T>> corr_me_;
    std::vector<int_dipole> insertion_terms_;
    std::vector<ab_term<T>> ab_terms_;
    std::vector<T> terms2_;
    T alphas_power_;
};

template <class T, class M, class S, class C, class R, class P, class U, class D>
using int_dipoles_integrand_t = int_dipoles_integrand<T, std::decay_t<M>, std::decay_t<S>,
    std::decay_t<C>, std::decay_t<R>, std::decay_t<P>, std::decay_t<U>, std::decay_t<D>>;

template <class T, class M, class S, class C, class R, class P, class U, class D>
inline std::unique_ptr<ps_integrand<T>> make_int_dipoles_integrand(
    M&& matrix_elements,
    S&& subtraction,
    C&& cuts,
    R&& recombiner,
    P&& pdfs,
    U&& scale_setter,
    D&& distributions,
    initial_state_set set,
    T hbarc2,
    finite_parts parts
) {
    return std::make_unique<int_dipoles_integrand_t<T, M, S, C, R, P, U, D>>(
        std::forward<M>(matrix_elements),
        std::forward<S>(subtraction),
        std::forward<C>(cuts),
        std::forward<R>(recombiner),
        std::forward<P>(pdfs),
        std::forward<U>(scale_setter),
        std::forward<D>(distributions),
        set,
        hbarc2,
        parts
    );
}

}

#endif
