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

		recombined_ps_.reserve(4 * (fs + 2));
		recombined_states_.reserve(fs);
		pdf_results_.reserve((pdfs_.count() == 1) ? 0 : pdfs_.count());

		if (!dynamic_scales_)
		{
			set_scales(std::vector<T>());
			results_.reserve(scales_.size());
		}

		pdfs_.register_partons(partons_in_initial_state_set(set));
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

		if (dynamic_scales_)
		{
			set_scales(recombined_ps_);
			results_.reserve(scales_.size());
		}

		T const factor = T(0.5) * hbarc2_ / info.energy_squared();

		borns_.resize(scales_.size());

		pdfs_.eval(info.x1(), scales_, scale_pdf_x1_, pdf_pdf_x1_);
		pdfs_.eval(info.x2(), scales_, scale_pdf_x2_, pdf_pdf_x2_);

		assert( scale_pdf_x1_.size() == scales_.size() );
		assert( scale_pdf_x2_.size() == scales_.size() );
		assert( (pdfs_.count() == 1) || (pdf_pdf_x1_.size() == pdfs_.count()) );
		assert( (pdfs_.count() == 1) || (pdf_pdf_x2_.size() == pdfs_.count()) );

		matrix_elements_.borns(phase_space, set_, scales_, borns_);

		assert( borns_.size() == scales_.size() );

		convolute_mes_with_pdfs(
			results_,
			pdf_results_,
			scale_pdf_x1_,
			scale_pdf_x2_,
			pdf_pdf_x1_,
			pdf_pdf_x2_,
			borns_,
			set_,
			factors_,
			factor,
			cut_result
		);

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

	std::vector<T> recombined_ps_;
	std::vector<final_state> final_states_;
	std::vector<recombined_state> recombined_states_;
	std::vector<parton_array<T>> scale_pdf_x1_;
	std::vector<parton_array<T>> scale_pdf_x2_;
	std::vector<parton_array<T>> pdf_pdf_x1_;
	std::vector<parton_array<T>> pdf_pdf_x2_;
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
