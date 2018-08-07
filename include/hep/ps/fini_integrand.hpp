#ifndef HEP_PS_FINI_INTEGRAND_HPP
#define HEP_PS_FINI_INTEGRAND_HPP

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

#include "hep/mc/projector.hpp"

#include "hep/ps/ab_terms.hpp"
#include "hep/ps/convolute.hpp"
#include "hep/ps/finite_parts.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/pdg_functions.hpp"
#include "hep/ps/ps_integrand.hpp"
#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scales.hpp"

#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace hep
{

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
class fini_integrand : public ps_integrand<T>
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
	fini_integrand(
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
		, insertion_terms_(matrix_elements_.insertion_terms())
		, alphas_power_(matrix_elements_.alphas_power())
	{
		std::size_t const fs = final_states_.size();

		recombined_ps_.reserve(4 * (fs + 2));
		recombined_states_.reserve(fs);
		pdf_results_.reserve((pdfs_.count() == 1) ? 0 : pdfs_.count());

		if (!scale_setter_.dynamic())
		{
			set_scales(std::vector<T>());
			results_.reserve(scales_.size());
			corr_me_.resize(insertion_terms_.size());
		}

		pdfs_.register_partons(partons_in_initial_state_set(set));

		for (auto const state : set)
		{
			me_set_.add(state);

			auto const one = state_parton_one(state);
			auto const two = state_parton_one(state);

			if (parton_type_of(one) == parton_type::quark)
			{
				me_set_.add(partons_to_initial_state(parton::photon, two));
				me_set_.add(partons_to_initial_state(parton::gluon, two));
			}

			if (parton_type_of(two) == parton_type::quark)
			{
				me_set_.add(partons_to_initial_state(one, parton::photon));
				me_set_.add(partons_to_initial_state(one, parton::gluon));
			}

			// FIXME: add missing cases
		}
	}

	T eval(
		std::vector<T> const& phase_space,
		luminosity_info<T> const& info,
		hep::projector<T>& projector
	) override {
		recombiner_.recombine(
			phase_space,
			final_states_,
			recombined_ps_,
			recombined_states_
		);

		auto const cut_result = cuts_.cut(
			recombined_ps_,
			info.rapidity_shift(),
			recombined_states_
		);

		if (cut_result.neg_cutted() && cut_result.pos_cutted())
		{
			return T();
		}

		if (scale_setter_.dynamic())
		{
			set_scales(recombined_ps_);
			results_.reserve(scales_.size());
			corr_me_.resize(insertion_terms_.size());
		}

		results_.clear();
		pdf_results_.clear();
		results_.resize(scales_.size());
		pdf_results_.resize((pdfs_.count() == 1) ? 0 : pdfs_.count());

		T const x = phase_space.back();
		T const x1p = info.x1() / (info.x1() * (T(1.0) - x) + x);
		T const x2p = info.x2() / (info.x2() * (T(1.0) - x) + x);

		pdfs_.eval(info.x1(), scales_, pdfsa1_, pdf_pdfsa1_);
		pdfs_.eval(info.x2(), scales_, pdfsa2_, pdf_pdfsa2_);
		pdfs_.eval(x1p, scales_, pdfsb1_, pdf_pdfsb1_);
		pdfs_.eval(x2p, scales_, pdfsb2_, pdf_pdfsb2_);

		assert( pdfsa1_.size() == scales_.size() );
		assert( pdfsa2_.size() == scales_.size() );
		assert( pdfsb1_.size() == scales_.size() );
		assert( pdfsb2_.size() == scales_.size() );
		assert( (pdfs_.count() == 1) || (pdf_pdfsa1_.size() == pdfs_.count()) );
		assert( (pdfs_.count() == 1) || (pdf_pdfsa2_.size() == pdfs_.count()) );
		assert( (pdfs_.count() == 1) || (pdf_pdfsb1_.size() == pdfs_.count()) );
		assert( (pdfs_.count() == 1) || (pdf_pdfsb2_.size() == pdfs_.count()) );

		for (auto& me : corr_me_)
		{
			me.clear();
		}

		matrix_elements_.correlated_me(phase_space, me_set_, corr_me_);
		auto const factor = T(0.5) * hbarc2_ / info.energy_squared();

		T const eta[] = { info.x1(), info.x2() };
		T const xprime[] = {
			eta[0] * (T(1.0) - x) + x,
			eta[1] * (T(1.0) - x) + x
		};

		// loop over all Born, FI, IF, II, and FF
		for (std::size_t index = 0; index != insertion_terms_.size(); ++index)
		{
			auto const me = corr_me_.at(index);
			auto const& term = insertion_terms_.at(index);

			if ((parts_ != finite_parts::insertion_term2) &&
				(term.type() != insertion_term_type::final_final))
			{
				std::size_t const i = term.initial_particle();

				subtraction_.insertion_terms(
					term,
					scales_,
					phase_space,
					xprime[1 - i],
					eta[1 - i],
					ab_neg_
				);

				assert( ab_neg_.size() == scales_.size() );

				subtraction_.insertion_terms(
					term,
					scales_,
					phase_space,
					xprime[i],
					eta[i],
					ab_pos_
				);

				assert( ab_neg_.size() == scales_.size() );

				for (std::size_t j = 0; j != scales_.size(); ++j)
				{
					auto const pdf_neg = effective_pdf(
						ab_neg_.at(j),
						xprime[1 - i],
						(i == 0) ? info.x2() : info.x1(),
						(i == 0) ? pdfsa2_.at(j) : pdfsa1_.at(j),
						(i == 0) ? pdfsb2_.at(j) : pdfsb1_.at(j)
					);

					auto const pdf_pos = effective_pdf(
						ab_pos_.at(j),
						xprime[i],
						(i == 0) ? info.x1() : info.x2(),
						(i == 0) ? pdfsa1_.at(j) : pdfsa2_.at(j),
						(i == 0) ? pdfsb1_.at(j) : pdfsb2_.at(j)
					);

					results_.at(j) += convolute(
						(i == 0) ? pdf_neg       : pdfsa2_.at(j),
						(i == 0) ? pdfsa1_.at(j) : pdf_neg,
						(i == 0) ? pdf_pos       : pdfsa1_.at(j),
						(i == 0) ? pdfsa2_.at(j) : pdf_pos,
						me,
						me_set_,
						factors_.at(j) * factor,
						cut_result
					);
				}

				for (std::size_t j = 0; j != pdf_pdfsa1_.size(); ++j)
				{
					auto const pdf_neg = effective_pdf(
						ab_neg_.front(),
						xprime[1 - i],
						(i == 0) ? info.x2() : info.x1(),
						(i == 0) ? pdf_pdfsa2_.at(j) : pdf_pdfsa1_.at(j),
						(i == 0) ? pdf_pdfsb2_.at(j) : pdf_pdfsb1_.at(j)
					);

					auto const pdf_pos = effective_pdf(
						ab_pos_.front(),
						xprime[i],
						(i == 0) ? info.x1() : info.x2(),
						(i == 0) ? pdf_pdfsa1_.at(j) : pdf_pdfsa2_.at(j),
						(i == 0) ? pdf_pdfsb1_.at(j) : pdf_pdfsb2_.at(j)
					);

					pdf_results_.at(j) += convolute(
						(i == 0) ? pdf_neg             : pdf_pdfsa2_.at(j),
						(i == 0) ? pdf_pdfsa1_.at(j) : pdf_neg,
						(i == 0) ? pdf_pos             : pdf_pdfsa1_.at(j),
						(i == 0) ? pdf_pdfsa2_.at(j) : pdf_pos,
						me,
						me_set_,
						factor,
						cut_result
					);
				}
			}

			if ((parts_ != finite_parts::insertion_term) &&
				(term.type() != insertion_term_type::born))
			{
				subtraction_.insertion_terms2(
					term,
					scales_,
					phase_space,
					terms2_
				);

				assert( terms2_.size() == scales_.size() );

				for (std::size_t i = 0; i != pdfsa1_.size(); ++i)
				{
					results_.at(i) += convolute(
						pdfsa1_.at(i),
						pdfsa2_.at(i),
						me,
						me_set_,
						factors_.at(i) * factor * terms2_.at(i),
						cut_result
					);
				}

				for (std::size_t i = 0; i != pdf_pdfsa1_.size(); ++i)
				{
					pdf_results_.at(i) += convolute(
						pdf_pdfsa1_.at(i),
						pdf_pdfsa2_.at(i),
						me,
						me_set_,
						factor * terms2_.front(),
						cut_result
					);
				}
			}
		}

		distributions_(
			recombined_ps_,
			info.rapidity_shift(),
			cut_result,
			results_,
			pdf_results_,
			projector
		);

		return results_.front().neg + results_.front().pos;
	}

protected:
	parton_array<T> effective_pdf(
		ab_terms<T> const& ab,
		T xprime,
		T eta,
		parton_array<T> const& pdfa,
		parton_array<T> const& pdfb
	) {
		parton_array<T> pdf;

		for (auto const a : parton_list())
		{
			for (auto const ap : parton_list())
			{
				// quark/anti-quark diagonality
				if ((a != ap) && ((a != parton::gluon) && (ap != parton::gluon)
					&& (a != parton::photon) && (ap != parton::photon)))
				{
					continue;
				}

				auto const at = parton_type_of(a);
				auto const apt = parton_type_of(ap);

				T q2 = T(1.0);

				if (((apt == parton_type::anti_quark) ||
					(apt == parton_type::quark)) &&
					(at == parton_type::photon_))
				{
					T const charge = T(pdg_id_to_charge_times_three(
						parton_to_pdg_id(ap))) / T(3.0);
					q2 = charge * charge;
				}

				pdf[a] += pdfb[ap] * q2 * ab.a[apt][at] * (T(1.0)-eta) / xprime;
				pdf[a] += pdfa[ap] * q2 * ab.b[apt][at];
			}
		}

		return pdf;
	}

	void set_scales(std::vector<T> const& phase_space)
	{
		using std::pow;

		scales_.clear();
		factors_.clear();
		scale_setter_(phase_space, scales_);
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

	std::vector<T> recombined_ps_;
	std::vector<final_state> final_states_;
	std::vector<recombined_state> recombined_states_;
	std::vector<parton_array<T>> pdfsa1_;
	std::vector<parton_array<T>> pdfsa2_;
	std::vector<parton_array<T>> pdfsb1_;
	std::vector<parton_array<T>> pdfsb2_;
	std::vector<parton_array<T>> pdf_pdfsa1_;
	std::vector<parton_array<T>> pdf_pdfsa2_;
	std::vector<parton_array<T>> pdf_pdfsb1_;
	std::vector<parton_array<T>> pdf_pdfsb2_;
	std::vector<neg_pos_results<T>> pdf_results_;
	std::vector<neg_pos_results<T>> results_;
	std::vector<scales<T>> scales_;
	std::vector<T> factors_;
	std::vector<initial_state_map<T>> corr_me_;
	std::vector<insertion_term> insertion_terms_;
	std::vector<ab_terms<T>> ab_neg_;
	std::vector<ab_terms<T>> ab_pos_;
	std::vector<T> terms2_;
	T alphas_power_;
};

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
using fini_integrand_t = fini_integrand<T, std::decay_t<M>, std::decay_t<S>,
	std::decay_t<C>, std::decay_t<R>, std::decay_t<P>, std::decay_t<U>,
	std::decay_t<D>>;

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
inline std::unique_ptr<ps_integrand<T>> make_fini_integrand(
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
	return std::make_unique<fini_integrand_t<T, M, S, C, R, P, U, D>>(
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
