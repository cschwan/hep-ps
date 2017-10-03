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
	{
		results_.reserve(pdfs_.count());
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

		pdfs_.eval(info.x1(), scales.factorization(), pdfsx1_);
		pdfs_.eval(info.x2(), scales.factorization(), pdfsx2_);

		auto const borns = matrix_elements_.borns(phase_space, set_);
		auto const factor = T(0.5) * hbarc2_ / info.energy_squared();

		std::size_t const size = pdfsx1_.size();
		results_.clear();

		for (std::size_t pdf = 0; pdf != size; ++pdf)
		{
			results_.push_back(convolute(
				pdfsx1_.at(pdf),
				pdfsx2_.at(pdf),
				borns,
				set_,
				factor,
				cut_result
			));
		}

		distributions_(phase_space, cut_result, results_, rapidity_shift,
			event_type::born_like_n, projector);

		return results_.front().neg + results_.front().pos;
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

	T old_renormalization_scale_;
	std::vector<parton_array<T>> pdfsx1_;
	std::vector<parton_array<T>> pdfsx2_;
	std::vector<neg_pos_results<T>> results_;
};

template <class T, class M, class C, class R, class P, class S, class D>
using born_integrand_t = born_integrand<T,
	typename std::decay<M>::type, typename std::decay<C>::type,
	typename std::decay<R>::type, typename std::decay<P>::type,
	typename std::decay<S>::type, typename std::decay<D>::type>;

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
	return std::unique_ptr<born_integrand_t<T, M, C, R, P, S, D>>(
		new born_integrand_t<T, M, C, R, P, S, D>(
			std::forward<M>(matrix_elements),
			std::forward<C>(cuts),
			std::forward<R>(recombiner),
			std::forward<P>(pdfs),
			std::forward<S>(scale_setter),
			std::forward<D>(distributions),
			set,
			hbarc2
	));
}

}

#endif
