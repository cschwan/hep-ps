#ifndef HEP_PS_FINI_INTEGRAND_HPP
#define HEP_PS_FINI_INTEGRAND_HPP

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

#include "hep/mc/projector.hpp"

#include "hep/ps/abc_terms.hpp"
#include "hep/ps/convolute.hpp"
#include "hep/ps/event_type.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/ps_integrand.hpp"
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
		bool insertion2
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
		, insertion2_(insertion2)
		, pdf_pdfsa1_(pdfs_.count())
		, pdf_pdfsa2_(pdfs_.count())
		, pdf_pdfsb1_(pdfs_.count())
		, pdf_pdfsb2_(pdfs_.count())
		, pdf_results_(pdfs_.count())
		, alphas_power_(matrix_elements_.alphas_power())
	{
		using std::begin;
		using std::end;

		auto const terms = matrix_elements_.insertion_terms();
		insertion_terms_.assign(begin(terms), end(terms));

		if (!scale_setter_.dynamic())
		{
			set_scales(std::vector<T>());
		}
	}

	T eval(
		std::vector<T> const& phase_space,
		luminosity_info<T> const& info,
		hep::projector<T>& projector
	) override {
		std::vector<T> aux_phase_space(phase_space.size());

		auto const recombined = recombiner_.recombine(
			phase_space,
			aux_phase_space,
			matrix_elements_.born_recombination_candidates(),
			0
		);

		if (recombined != 0)
		{
			return T();
		}

		T const rapidity_shift = info.rapidity_shift();
		auto const cut_result = cuts_.cut(phase_space, rapidity_shift,
			event_type::born_like_n);

		if (cut_result.neg_cutted() && cut_result.pos_cutted())
		{
			return T();
		}

		if (scale_setter_.dynamic())
		{
			set_scales(phase_space);
		}

		pdfsa1_.clear();
		pdfsa1_.resize(scales_.size());
		pdfsa2_.clear();
		pdfsa2_.resize(scales_.size());
		pdfsb1_.clear();
		pdfsb1_.resize(scales_.size());
		pdfsb2_.clear();
		pdfsb2_.resize(scales_.size());

		T const x = phase_space.back();
		T const x1p = info.x1() / (info.x1() * (T(1.0) - x) + x);
		T const x2p = info.x2() / (info.x2() * (T(1.0) - x) + x);

		pdfs_.eval(info.x1(), scales_, pdfsa1_);
		pdfs_.eval(info.x2(), scales_, pdfsa2_);
		pdfs_.eval(x1p, scales_, pdfsb1_);
		pdfs_.eval(x2p, scales_, pdfsb2_);

		results_.clear();
		results_.resize(scales_.size());

		if (pdfs_.count() > 1)
		{
			pdf_pdfsa1_.clear();
			pdf_pdfsa1_.resize(pdfs_.count());
			pdf_pdfsa2_.clear();
			pdf_pdfsa2_.resize(pdfs_.count());
			pdf_pdfsb1_.clear();
			pdf_pdfsb1_.resize(pdfs_.count());
			pdf_pdfsb2_.clear();
			pdf_pdfsb2_.resize(pdfs_.count());

			T const muf = scales_.front().factorization();

			pdfs_.eval(info.x1(), muf, pdf_pdfsa1_);
			pdfs_.eval(info.x2(), muf, pdf_pdfsa2_);
			pdfs_.eval(x1p, muf, pdf_pdfsb1_);
			pdfs_.eval(x2p, muf, pdf_pdfsb2_);

			pdf_results_.clear();
			pdf_results_.resize(pdfs_.count());
		}

		auto const& corr_me = matrix_elements_.correlated_me(phase_space, set_);
		auto const factor = T(0.5) * hbarc2_ / info.energy_squared();

		// loop over all Born, FI, IF, and II
		for (std::size_t index = 0; index != insertion_terms_.size(); ++index)
		{
			auto const me = corr_me.at(index);
			auto const& term = insertion_terms_.at(index);

			// loop over both initial state partons
			for (auto const i : { 0u, 1u })
			{
				switch (term.type())
				{
				case insertion_term_type::final_initial:
					if (term.spectator() != i)
					{
						continue;
					}

					break;

				case insertion_term_type::initial_final:
				case insertion_term_type::initial_initial:
					if (term.emitter() != i)
					{
						continue;
					}

					break;

				default:
					break;
				}

				T const eta_neg = (i == 0) ? info.x2() : info.x1();
				T const xprime_neg = eta_neg * (T(1.0) - x) + x;

				abc_neg_.clear();
				subtraction_.insertion_terms(
					term,
					scales_,
					phase_space,
					xprime_neg,
					eta_neg,
					abc_neg_
				);

				T const eta_pos = (i == 0) ? info.x1() : info.x2();
				T const xprime_pos = eta_pos * (T(1.0) - x) + x;

				abc_pos_.clear();
				subtraction_.insertion_terms(
					term,
					scales_,
					phase_space,
					xprime_pos,
					eta_pos,
					abc_pos_
				);

				for (std::size_t j = 0; j != scales_.size(); ++j)
				{
					auto const pdf_neg = effective_pdf(
						abc_neg_.at(j),
						xprime_neg,
						(i == 0) ? info.x2() : info.x1(),
						(i == 0) ? pdfsa2_.at(j) : pdfsa1_.at(j),
						(i == 0) ? pdfsb2_.at(j) : pdfsb1_.at(j)
					);

					auto const pdf_pos = effective_pdf(
						abc_pos_.at(j),
						xprime_pos,
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
						set_,
						factors_.at(j) * factor,
						cut_result
					);
				}

				if (pdfs_.count() == 1)
				{
					continue;
				}

				for (std::size_t pdf = 0; pdf != pdfs_.count(); ++pdf)
				{
					auto const pdf_neg = effective_pdf(
						abc_neg_.front(),
						xprime_neg,
						(i == 0) ? info.x2() : info.x1(),
						(i == 0) ? pdf_pdfsa2_.at(pdf) : pdf_pdfsa1_.at(pdf),
						(i == 0) ? pdf_pdfsb2_.at(pdf) : pdf_pdfsb1_.at(pdf)
					);

					auto const pdf_pos = effective_pdf(
						abc_pos_.front(),
						xprime_pos,
						(i == 0) ? info.x1() : info.x2(),
						(i == 0) ? pdf_pdfsa1_.at(pdf) : pdf_pdfsa2_.at(pdf),
						(i == 0) ? pdf_pdfsb1_.at(pdf) : pdf_pdfsb2_.at(pdf)
					);

					pdf_results_.at(pdf) += convolute(
						(i == 0) ? pdf_neg             : pdf_pdfsa2_.at(pdf),
						(i == 0) ? pdf_pdfsa1_.at(pdf) : pdf_neg,
						(i == 0) ? pdf_pos             : pdf_pdfsa1_.at(pdf),
						(i == 0) ? pdf_pdfsa2_.at(pdf) : pdf_pos,
						me,
						set_,
						factor,
						cut_result
					);
				}
			}

			if (insertion2_)
			{
				subtraction_.insertion_terms2(
					term,
					scales_,
					phase_space,
					terms2_
				);

				for (std::size_t i = 0; i != scales_.size(); ++i)
				{
					results_.at(i) += convolute(
						pdfsa1_.front(),
						pdfsa2_.front(),
						me,
						set_,
						factors_.at(i) * factor * terms2_.at(i),
						cut_result
					);
				}

				if (pdfs_.count() == 1)
				{
					continue;
				}

				for (std::size_t pdf = 0; pdf != pdfs_.count(); ++pdf)
				{
					pdf_results_.at(pdf) += convolute(
						pdf_pdfsa1_.at(pdf),
						pdf_pdfsa2_.at(pdf),
						me,
						set_,
						factor * terms2_.front(),
						cut_result
					);
				}
			}
		}

		distributions_(
			phase_space,
			cut_result,
			results_,
			pdf_results_,
			rapidity_shift,
			event_type::born_like_n,
			projector
		);

		return results_.front().neg + results_.front().pos;
	}

protected:
	parton_array<T> effective_pdf(
		abc_terms<T> const& abc,
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
				if ((a != parton::gluon) && (ap != parton::gluon) &&
					(a != ap))
				{
					continue;
				}

				auto const at = parton_type_of(a);
				auto const apt = parton_type_of(ap);

				pdf[a] += pdfb[ap] * abc.a[apt][at] * (T(1.0) - eta) / xprime;
				pdf[a] -= pdfa[ap] * abc.b[apt][at] * (T(1.0) - eta);
				pdf[a] += pdfa[ap] * abc.c[apt][at];
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
	T hbarc2_;

	bool insertion2_;
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
	std::vector<insertion_term> insertion_terms_;
	std::vector<abc_terms<T>> abc_neg_;
	std::vector<abc_terms<T>> abc_pos_;
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
	bool insertion2 = false
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
		insertion2
	);
}

}

#endif
