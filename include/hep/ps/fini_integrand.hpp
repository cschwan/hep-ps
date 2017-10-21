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
		, pdfsa1_(pdfs_.count())
		, pdfsa2_(pdfs_.count())
		, pdfsb1_(pdfs_.count())
		, pdfsb2_(pdfs_.count())
		, results_(pdfs_.count())
	{
		using std::begin;
		using std::end;

		auto const terms = matrix_elements_.insertion_terms();
		insertion_terms_.assign(begin(terms), end(terms));
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

		auto const scales = scale_setter_(phase_space);

		// only set renormalization scale if it changed
		if (scales.renormalization() != old_renormalization_scale_)
		{
			old_renormalization_scale_ = scales.renormalization();

			T const alphas = pdfs_.alphas(scales.renormalization());
			matrix_elements_.parameters(scales.renormalization(), alphas);
		}

		T const x = phase_space.back();
		T const muf = scales.factorization();

		pdfs_.eval(info.x1(), muf, pdfsa1_);
		pdfs_.eval(info.x2(), muf, pdfsa2_);
		pdfs_.eval(info.x1() / (info.x1() * (T(1.0) - x) + x), muf, pdfsb1_);
		pdfs_.eval(info.x2() / (info.x2() * (T(1.0) - x) + x), muf, pdfsb2_);

		auto const& corr_me = matrix_elements_.correlated_me(phase_space, set_);
		auto const factor = T(0.5) * hbarc2_ / info.energy_squared();

		auto effective_pdf = [&](
			insertion_term const& term,
			T eta,
			parton_array<T> const& pdfa,
			parton_array<T> const& pdfb
		) {
			T const ome = (T(1.0) - eta);
			T const xprime = eta + ome * x;
			auto const& abc = subtraction_.insertion_terms(
				term,
				scales,
				phase_space,
				xprime,
				eta
			);

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

					pdf[a] += pdfb[ap] * abc.a[apt][at] * ome / xprime;
					pdf[a] -= pdfa[ap] * abc.b[apt][at] * ome;
					pdf[a] += pdfa[ap] * abc.c[apt][at];
				}
			}

			return pdf;
		};

		std::size_t const size = pdfsa1_.size();
		results_.clear();
		results_.resize(size);

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

				for (std::size_t pdf = 0; pdf != size; ++pdf)
				{
					auto const pdf_neg = effective_pdf(
						term,
						(i == 0) ? info.x2() : info.x1(),
						(i == 0) ? pdfsa2_.at(pdf) : pdfsa1_.at(pdf),
						(i == 0) ? pdfsb2_.at(pdf) : pdfsb1_.at(pdf)
					);

					auto const pdf_pos = effective_pdf(
						term,
						(i == 0) ? info.x1() : info.x2(),
						(i == 0) ? pdfsa1_.at(pdf) : pdfsa2_.at(pdf),
						(i == 0) ? pdfsb1_.at(pdf) : pdfsb2_.at(pdf)
					);

					if (i == 0)
					{
						results_.at(pdf) += convolute(
							pdf_neg,
							pdfsa1_.at(pdf),
							pdf_pos,
							pdfsa2_.at(pdf),
							me,
							set_,
							factor,
							cut_result
						);
					}
					else
					{
						results_.at(pdf) += convolute(
							pdfsa2_.at(pdf),
							pdf_neg,
							pdfsa1_.at(pdf),
							pdf_pos,
							me,
							set_,
							factor,
							cut_result
						);
					}
				}
			}

			if (insertion2_)
			{
				T const ins = factor * subtraction_.insertion_terms2(
					term,
					scales,
					phase_space
				);

				for (std::size_t pdf = 0; pdf != size; ++pdf)
				{
					results_.at(pdf) += convolute(pdfsa1_.at(pdf),
						pdfsa2_.at(pdf), me, set_, ins, cut_result);
				}
			}
		}

		distributions_(phase_space, cut_result, results_, rapidity_shift,
			event_type::born_like_n, projector);

		return results_.front().neg + results_.front().pos;
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

	T old_renormalization_scale_;
	bool insertion2_;
	std::vector<parton_array<T>> pdfsa1_;
	std::vector<parton_array<T>> pdfsa2_;
	std::vector<parton_array<T>> pdfsb1_;
	std::vector<parton_array<T>> pdfsb2_;
	std::vector<neg_pos_results<T>> results_;
	std::vector<insertion_term> insertion_terms_;
};

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
using fini_integrand_t = fini_integrand<T,
	typename std::decay<M>::type, typename std::decay<S>::type,
	typename std::decay<C>::type, typename std::decay<R>::type,
	typename std::decay<P>::type, typename std::decay<U>::type,
	typename std::decay<D>::type>;

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
	return std::unique_ptr<fini_integrand_t<T, M, S, C, R, P, U, D>>(
		new fini_integrand_t<T, M, S, C, R, P, U, D>(
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
	));
}

}

#endif
