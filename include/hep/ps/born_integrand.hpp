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
		, pdf_uncertainity_results_(pdfs_.count() - 1)
		, alphas_power_(matrix_elements_.alphas_power())
		, dynamic_scales_{scale_setter_.dynamic() && (alphas_power_ > T{})}
	{
		pdf_uncertainity_results_.reserve(pdfs_.count());

		if (!dynamic_scales_)
		{
			auto const central_scale = scale_setter_(std::vector<T>(), scales_);
			alphas_ = pdfs_.eval_alphas(central_scale.renormalization());
			matrix_elements_.alphas(alphas_);

			for (auto const scale : scales_)
			{
				using std::pow;

				T const alphas = pdfs_.eval_alphas(scale.renormalization());
				static_factors_.push_back(pow(alphas / alphas_, alphas_power_));
			}
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

		auto const central_scale = scale_setter_(phase_space, scales_);

		if (dynamic_scales_)
		{
			alphas_ = pdfs_.eval_alphas(central_scale.renormalization());
			matrix_elements_.alphas(alphas_);
		}

		auto const factor = T(0.5) * hbarc2_ / info.energy_squared();

		// for the central scale calculate scale-dependent and scale-independent
		// parts separately
		auto borns = matrix_elements_.borns(phase_space, set_);
		auto const scale_dep_me = matrix_elements_.borns(phase_space, set_,
			central_scale.renormalization());

		scale_uncertainty_results_.clear();

		// evaluate PDFs for the non-central scales
		for (std::size_t i = 0; i != scales_.size(); ++i)
		{
			T const muf = scales_.at(i).factorization();
			T const mur = scales_.at(i).renormalization();

			// TODO: avoid evaluating the same PDFs
			auto const pdfx1 = pdfs_.eval(info.x1(), muf);
			auto const pdfx2 = pdfs_.eval(info.x2(), muf);

			// rescaling factor from the strong coupling
			T rescaling = T(1.0);

			if (dynamic_scales_)
			{
				using std::pow;

				// TODO: avoid evaluating the same alphas
				T const alphas = pdfs_.eval_alphas(mur);
				rescaling = pow(alphas / alphas_, alphas_power_);
			}
			else
			{
				rescaling = static_factors_.at(i);
			}

			// TODO: avoid calculating the same matrix elements
			auto new_borns = matrix_elements_.borns(phase_space, set_, mur);
			new_borns += borns;

			scale_uncertainty_results_.push_back(convolute(
				pdfx1,
				pdfx2,
				new_borns,
				set_,
				factor * rescaling,
				cut_result
			));
		}

		// now we only need the result for the central scale
		borns += scale_dep_me;

		// evaluate all PDFs for the central scale
		pdfs_.eval(info.x1(), central_scale.factorization(), pdfsx1_);
		pdfs_.eval(info.x2(), central_scale.factorization(), pdfsx2_);

		pdf_uncertainity_results_.clear();

		for (std::size_t pdf = 1; pdf < pdfsx1_.size(); ++pdf)
		{
			pdf_uncertainity_results_.push_back(convolute(
				pdfsx1_.at(pdf),
				pdfsx2_.at(pdf),
				borns,
				set_,
				factor,
				cut_result
			));
		}

		// calculate the central result
		auto const result = convolute(
			pdfsx1_.front(),
			pdfsx2_.front(),
			borns,
			set_,
			factor,
			cut_result
		);

		distributions_(
			phase_space,
			cut_result,
			result,
			scale_uncertainty_results_,
			pdf_uncertainity_results_,
			rapidity_shift,
			event_type::born_like_n,
			projector
		);

		return result.neg + result.pos;
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
	std::vector<neg_pos_results<T>> pdf_uncertainity_results_;
	std::vector<neg_pos_results<T>> scale_uncertainty_results_;
	std::vector<scales<T>> scales_;
	std::vector<T> static_factors_;
	T alphas_power_;
	T alphas_;
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
