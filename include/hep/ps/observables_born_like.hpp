#ifndef HEP_PS_OBSERVABLES_BORN_LIKE_HPP
#define HEP_PS_OBSERVABLES_BORN_LIKE_HPP

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

#include "hep/ps/event_type.hpp"
#include "hep/ps/fold.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/particle_type.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <utility>
#include <type_traits>
#include <vector>

namespace hep
{

template <class T, class M, class C, class R, class P, class S, class D>
class observables_born_like
{
public:
	template <
		typename MatrixElements,
		typename Cuts,
		typename Recombiner,
		typename Pdf,
		typename ScaleSetter,
		typename Distributions>
	observables_born_like(
		MatrixElements&& matrix_elements,
		Cuts&& cuts,
		Recombiner&& recombiner,
		Pdf&& pdf,
		ScaleSetter&& scale_setter,
		Distributions&& distributions,
		T hbarc2
	)
		: matrix_elements_(std::forward<MatrixElements>(matrix_elements))
		, cuts_(std::forward<Cuts>(cuts))
		, recombiner_(std::forward<Recombiner>(recombiner))
		, pdf_(std::forward<Pdf>(pdf))
		, scale_setter_(std::forward<ScaleSetter>(scale_setter))
		, distributions_(std::forward<Distributions>(distributions))
		, hbarc2_(hbarc2)
	{
	}

	T operator()(
		std::vector<T> const& phase_space,
		luminosity_info<T> const& info,
		initial_state_set set
	) {
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
			matrix_elements_.scale(scales.renormalization(), pdf_);
			old_renormalization_scale_ = scales.renormalization();
		}

		auto const borns = matrix_elements_.borns(phase_space, set);
		auto const pdfx1 = pdf_.pdf(info.x1(), scales.factorization());
		auto const pdfx2 = pdf_.pdf(info.x2(), scales.factorization());
		auto const factor = T(0.5) * hbarc2_ / info.energy_squared();
		auto const result = fold(pdfx1, pdfx2, borns, set, factor, cut_result);

		distributions_(phase_space, cut_result, result, rapidity_shift,
			event_type::born_like_n);

		return result.neg + result.pos;
	}

	D const& distributions() const
	{
		return distributions_;
	}

	D& distributions()
	{
		return distributions_;
	}

	M const& matrix_elements() const
	{
		return matrix_elements_;
	}

	M& matrix_elements()
	{
		return matrix_elements_;
	}

private:
	M matrix_elements_;
	C cuts_;
	R recombiner_;
	P pdf_;
	S scale_setter_;
	D distributions_;
	T hbarc2_;

	T old_renormalization_scale_;
};

template <class T, class M, class C, class R, class P, class S, class D>
using observables_born_like_type = observables_born_like<T,
	typename std::decay<M>::type, typename std::decay<C>::type,
	typename std::decay<R>::type, typename std::decay<P>::type,
	typename std::decay<S>::type, typename std::decay<D>::type>;

template <class T, class M, class C, class R, class P, class S, class D>
inline observables_born_like_type<T, M, C, R, P, S, D>
make_observables_born_like(
	M&& matrix_elements,
	C&& cuts,
	R&& recombiner,
	P&& pdf,
	S&& scale_setter,
	D&& distributions,
	T hbarc2
) {
	return observables_born_like_type<T, M, C, R, P, S, D>(
		std::forward<M>(matrix_elements),
		std::forward<C>(cuts),
		std::forward<R>(recombiner),
		std::forward<P>(pdf),
		std::forward<S>(scale_setter),
		std::forward<D>(distributions),
		hbarc2
	);
}

}

#endif
