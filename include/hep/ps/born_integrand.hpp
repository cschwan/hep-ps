#ifndef HEP_PS_BORN_INTEGRAND_HPP
#define HEP_PS_BORN_INTEGRAND_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2016-2017  Christopher Schwan
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
#include "hep/ps/event_type.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/ps_integrand.hpp"
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
		, pdfsx1_(pdfs_.count())
		, pdfsx2_(pdfs_.count())
		, pdf_results_(pdfs_.count() - 1)
		, alphas_power_(matrix_elements_.alphas_power())
		, dynamic_scales_{scale_setter_.dynamic()}
	{
		pdf_results_.reserve(pdfs_.count());

		if (!dynamic_scales_)
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

		if (dynamic_scales_)
		{
			set_scales(phase_space);
		}

		T const factor = T(0.5) * hbarc2_ / info.energy_squared();

		results_.clear();
		borns_.resize(scales_.size());
		pdfsx1_.resize(scales_.size());
		pdfsx2_.resize(scales_.size());

		pdfs_.eval(info.x1(), scales_, pdfsx1_);
		pdfs_.eval(info.x2(), scales_, pdfsx2_);

		matrix_elements_.borns(phase_space, set_, scales_, borns_);

		for (std::size_t i = 0; i != scales_.size(); ++i)
		{
			results_.push_back(convolute(
				pdfsx1_.at(i),
				pdfsx2_.at(i),
				borns_.at(i),
				set_,
				factors_.at(i) * factor,
				cut_result
			));
		}

		if (pdfs_.count() > 1)
		{
			pdf_results_.clear();
			pdfsx1_.resize(pdfs_.count());
			pdfsx2_.resize(pdfs_.count());

			// TODO: do not evaluate the central PDF, we don't need it

			// evaluate all PDFs for the central scale
			pdfs_.eval(info.x1(), scales_.front().factorization(), pdfsx1_);
			pdfs_.eval(info.x2(), scales_.front().factorization(), pdfsx2_);

			for (std::size_t pdf = 1; pdf != pdfsx1_.size(); ++pdf)
			{
				pdf_results_.push_back(convolute(
					pdfsx1_.at(pdf),
					pdfsx2_.at(pdf),
					borns_.at(0),
					set_,
					factor,
					cut_result
				));
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
	C cuts_;
	R recombiner_;
	P pdfs_;
	S scale_setter_;
	D distributions_;
	initial_state_set set_;
	T hbarc2_;

	std::vector<parton_array<T>> pdfsx1_;
	std::vector<parton_array<T>> pdfsx2_;
	std::vector<neg_pos_results<T>> pdf_results_;
	std::vector<neg_pos_results<T>> results_;
	std::vector<scales<T>> scales_;
	std::vector<T> factors_;
	std::vector<initial_state_array<T>> borns_;
	T alphas_power_;
	bool dynamic_scales_;
};

template <class T, class M, class C, class R, class P, class S, class D>
using born_integrand_t = born_integrand<T, std::decay_t<M>, std::decay_t<C>,
	std::decay_t<R>, std::decay_t<P>, std::decay_t<S>, std::decay_t<D>>;

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
