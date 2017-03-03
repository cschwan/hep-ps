#ifndef HEP_PS_OBSERVABLES_FINITE_INSERTION_HPP
#define HEP_PS_OBSERVABLES_FINITE_INSERTION_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
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
#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/trivial_distributions.hpp"

#include <cassert>
#include <cstddef>
#include <utility>
#include <type_traits>
#include <vector>

namespace hep
{

template <class T, class M, class S, class C, class R, class P, class U>
class observables_finite_insertion
{
public:
	template <
		typename MatrixElements,
		typename Subtraction,
		typename Cuts,
		typename Recombiner,
		typename Pdf,
		typename ScaleSetter>
	observables_finite_insertion(
		MatrixElements&& matrix_elements,
		Subtraction&& subtraction,
		Cuts&& cuts,
		Recombiner&& recombiner,
		Pdf&& pdf,
		ScaleSetter&& scale_setter,
		T hbarc2
	)
		: matrix_elements_(std::forward<MatrixElements>(matrix_elements))
		, subtraction_(std::forward<Subtraction>(subtraction))
		, cuts_(std::forward<Cuts>(cuts))
		, recombiner_(std::forward<Recombiner>(recombiner))
		, pdf_(std::forward<Pdf>(pdf))
		, scale_setter_(std::forward<ScaleSetter>(scale_setter))
		, hbarc2_(hbarc2)
	{
	}

	template <typename Distributions = trivial_distributions<T>>
	T operator()(
		std::vector<T> const& phase_space,
		luminosity_info<T> const& info,
		T x,
		initial_state_set set,
		Distributions&& distributions = trivial_distributions<T>()
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

		auto d_function = [&](
			std::size_t i,
			initial_state state,
			decltype (cut_result) cut,
			abc_terms<T> const& abc,
			initial_state_array<T> const& me
		) {
			auto const a = (i == 0)
				? state_parton_one(state)
				: state_parton_two(state);

			T d{};

			if (a == parton::gluon)
			{
				// TODO: NYI
				assert( false );
			}

			for (auto const ap : parton_list())
			{
				d += pdfb[i][ap] * abc.a[ap][a] * (T(1.0) - eta[i]) / xprime[i];
				d -= pdfa[i][ap] * abc.b[ap][a] * (T(1.0) - eta[i]);
				d += pdfa[i][ap] * abc.c[ap][a];
			}

			auto const b = (i == 0)
				? state_parton_two(state)
				: state_parton_one(state);
			auto const j = 1 - i;

			neg_pos_results<T> result;

			if (!cut.neg_cutted() && state_has_neg_shift(state))
			{
				result.neg = d * pdfa[j][b] * me[state];
			}

			if (!cut.pos_cutted() && state_has_pos_shift(state))
			{
				result.pos = d * pdfa[j][b] * me[state];
			}

			return result;
		};

		auto const borns = matrix_elements_.borns(phase_space, set);

		neg_pos_results<T> result;

		// loop over the both initial particles
		for (auto const i : { 0, 1 })
		{
			auto const abc = subtraction_.finite_born(xprime[i], eta[i]);

			// loop over all initial states
			for (auto const state : set)
			{
				result += d_function(i, state, cut_result, abc, borns);
			}
		}

		auto const& insertion_terms = matrix_elements_.insertion_terms();
		auto const mu2 = scales.factorization() * scales.factorization();
		auto const& correlated_me = matrix_elements_.correlated_me(phase_space,
			set);

		// loop over the two particle in the initial state
		for (auto const i : { 0, 1 })
		{
			// loop over all FI, IF, and II
			for (std::size_t index = 0; index != insertion_terms.size(); ++index)
			{
				auto const& term = insertion_terms.at(index);
				abc_terms<T> abc;

				switch (term.type())
				{
				case insertion_term_type::final_initial:
					abc = subtraction_.finite_final_initial(
						xprime[i],
						eta[i],
						term.emitter_type()
					);
					break;

				case insertion_term_type::initial_final:
					abc = subtraction_.finite_initial_final(
						xprime[i],
						eta[i],
						mu2,
						phase_space,
						term.emitter(),
						term.spectator()
					);
					break;

				case insertion_term_type::initial_initial:
					abc = subtraction_.finite_initial_initial(
						xprime[i],
						eta[i],
						mu2,
						phase_space,
						term.emitter(),
						term.spectator()
					);
					break;
				}

				auto const me = correlated_me.at(index);

				for (auto const state : set)
				{
					result += d_function(i, state, cut_result, abc, me);
				}
			}
		}

		distributions(phase_space, cut_result, result, rapidity_shift,
			event_type::born_like_n);

		return result.neg + result.pos;
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
	T hbarc2_;

	T old_renormalization_scale_;
};

template <class T, class M, class S, class C, class R, class P, class U>
using observables_finite_insertion_type = observables_finite_insertion<T,
	typename std::decay<M>::type, typename std::decay<S>::type,
	typename std::decay<C>::type, typename std::decay<R>::type,
	typename std::decay<P>::type, typename std::decay<U>::type>;

template <class T, class M, class S, class C, class R, class P, class U>
inline observables_finite_insertion_type<T, M, S, C, R, P, U>
make_observables_finite_insertion(
	M&& matrix_elements,
	S&& subtraction,
	C&& cuts,
	R&& recombiner,
	P&& pdf,
	U&& scale_setter,
	T hbarc2
) {
	return observables_finite_insertion_type<T, M, S, C, R, P, U>(
		std::forward<M>(matrix_elements),
		std::forward<S>(subtraction),
		std::forward<C>(cuts),
		std::forward<R>(recombiner),
		std::forward<P>(pdf),
		std::forward<U>(scale_setter),
		hbarc2
	);
}

}

#endif
