#ifndef HEP_PS_OBSERVABLES_FINI_HPP
#define HEP_PS_OBSERVABLES_FINI_HPP

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

#include "hep/ps/abc_terms.hpp"
#include "hep/ps/event_type.hpp"
#include "hep/ps/fold.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/observables.hpp"
#include "hep/ps/random_numbers.hpp"

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
class observables_fini : public observables<T>
{
public:
	template <
		typename MatrixElements,
		typename Subtraction,
		typename Cuts,
		typename Recombiner,
		typename Pdf,
		typename ScaleSetter,
		typename Distributions>
	observables_fini(
		MatrixElements&& matrix_elements,
		Subtraction&& subtraction,
		Cuts&& cuts,
		Recombiner&& recombiner,
		Pdf&& pdf,
		ScaleSetter&& scale_setter,
		Distributions&& distributions,
		T hbarc2,
		bool insertion2
	)
		: matrix_elements_(std::forward<MatrixElements>(matrix_elements))
		, subtraction_(std::forward<Subtraction>(subtraction))
		, cuts_(std::forward<Cuts>(cuts))
		, recombiner_(std::forward<Recombiner>(recombiner))
		, pdf_(std::forward<Pdf>(pdf))
		, scale_setter_(std::forward<ScaleSetter>(scale_setter))
		, distributions_(std::forward<Distributions>(distributions))
		, hbarc2_(hbarc2)
		, insertion2_(insertion2)
	{
	}

	T eval(
		std::vector<T> const& phase_space,
		luminosity_info<T> const& info,
		random_numbers<T>& extra_random_numbers,
		initial_state_set set
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
			matrix_elements_.scale(scales.renormalization(), pdf_);
			old_renormalization_scale_ = scales.renormalization();
		}

		T const x = extra_random_numbers.front();
		T const eta[] = {
			info.x1(),
			info.x2()
		};
		T const xprime[] = {
			eta[0] + (T(1.0) - eta[0]) * x,
			eta[1] + (T(1.0) - eta[1]) * x
		};
		parton_array<T> const pdfa[] = {
			pdf_.pdf(eta[0], scales.factorization()),
			pdf_.pdf(eta[1], scales.factorization())
		};
		parton_array<T> const pdfb[] = {
			pdf_.pdf(eta[0] / xprime[0], scales.factorization()),
			pdf_.pdf(eta[1] / xprime[1], scales.factorization())
		};

		auto const& corr_me = matrix_elements_.correlated_me(phase_space, set);
		auto const& insertion_terms = matrix_elements_.insertion_terms();
		auto const factor = T(0.5) * hbarc2_ / info.energy_squared();

		neg_pos_results<T> result;

		// loop over all Born, FI, IF, and II
		for (std::size_t index = 0; index != insertion_terms.size(); ++index)
		{
			auto const me = corr_me.at(index);

			// loop over both initial state partons
			for (auto const i : { 0, 1 })
			{
				auto const& abc = subtraction_.insertion_terms(
					insertion_terms.at(index),
					scales,
					phase_space,
					xprime[i],
					eta[i],
					i
				);

				parton_array<T> d;

				for (auto const a : parton_list())
				{
					for (auto const ap : parton_list())
					{
						d[a] += pdfb[i][ap] * abc.a[ap][a] * (T(1.0) - eta[i]) /
							xprime[i];
						d[a] -= pdfa[i][ap] * abc.b[ap][a] * (T(1.0) - eta[i]);
						d[a] += pdfa[i][ap] * abc.c[ap][a];
					}
				}

				if (i == 0)
				{
					result += fold(d, pdfa[1], me, set, factor, cut_result);
				}
				else
				{
					result += fold(pdfa[0], d, me, set, factor, cut_result);
				}
			}

			if (insertion2_)
			{
				T const ins = factor * subtraction_.insertion_terms2(
					insertion_terms.at(index),
					scales,
					phase_space
				);

				result += fold(pdfa[0], pdfa[1], me, set, ins, cut_result);
			}
		}

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
	S subtraction_;
	C cuts_;
	R recombiner_;
	P pdf_;
	U scale_setter_;
	D distributions_;
	T hbarc2_;

	T old_renormalization_scale_;
	bool insertion2_;
};

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
using observables_fini_t = observables_fini<T,
	typename std::decay<M>::type, typename std::decay<S>::type,
	typename std::decay<C>::type, typename std::decay<R>::type,
	typename std::decay<P>::type, typename std::decay<U>::type,
	typename std::decay<D>::type>;

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
inline std::unique_ptr<observables<T>> make_observables_fini(
	M&& matrix_elements,
	S&& subtraction,
	C&& cuts,
	R&& recombiner,
	P&& pdf,
	U&& scale_setter,
	D&& distributions,
	T hbarc2,
	bool insertion2 = false
) {
	return std::unique_ptr<observables_fini_t<T, M, S, C, R, P, U, D>>(
		new observables_fini_t<T, M, S, C, R, P, U, D>(
			std::forward<M>(matrix_elements),
			std::forward<S>(subtraction),
			std::forward<C>(cuts),
			std::forward<R>(recombiner),
			std::forward<P>(pdf),
			std::forward<U>(scale_setter),
			std::forward<D>(distributions),
			hbarc2,
			insertion2
	));
}

}

#endif
