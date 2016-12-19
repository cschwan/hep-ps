#ifndef HEP_PS_OBSERVABLES_BORN_LIKE_HPP
#define HEP_PS_OBSERVABLES_BORN_LIKE_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016  Christopher Schwan
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

#include "hep/ps/cut_result.hpp"
#include "hep/ps/initial_state_array.hpp"
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

namespace
{

// TODO: Make `constexpr` in C++14
inline bool cut_required(hep::initial_state state, hep::cut_result cut)
{
	return (hep::state_has_neg_shift(state) && cut.neg_cutted()) ||
		(hep::state_has_pos_shift(state) && cut.pos_cutted());
}

}

namespace hep
{

template <class T, class M, class C, class R, class L, class S>
class observables_born_like
{
public:
	template <
		typename MatrixElements,
		typename Cuts,
		typename Recombiner,
		typename Luminosities,
		typename ScaleSetter>
	observables_born_like(
		MatrixElements&& matrix_elements,
		Cuts&& cuts,
		Recombiner&& recombiner,
		Luminosities&& luminosities,
		ScaleSetter&& scale_setter,
		T conversion_constant
	)
		: matrix_elements_(std::forward<MatrixElements>(matrix_elements))
		, cuts_(std::forward<Cuts>(cuts))
		, recombiner_(std::forward<Recombiner>(recombiner))
		, luminosities_(std::forward<Luminosities>(luminosities))
		, scale_setter_(std::forward<ScaleSetter>(scale_setter))
		, conversion_constant_(conversion_constant)
	{
	}

	T operator()(
		std::vector<T> const& phase_space,
		luminosity_info<T> const& info,
		initial_state_set set
	) {
		// TODO: generate distributions

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
		auto const cut_result = cuts_.cut(phase_space, rapidity_shift, false);

		if (cut_result.neg_cutted() && cut_result.pos_cutted())
		{
			return T();
		}

		auto const scales = scale_setter_(phase_space);

		// only set renormalization scale if it changed
		if (scales.renormalization() != old_renormalization_scale_)
		{
			matrix_elements_.scale(scales.renormalization(), luminosities_);
			old_renormalization_scale_ = scales.renormalization();
		}

		auto borns = matrix_elements_.borns(phase_space, set);

		if (cut_result.neg_cutted() || cut_result.pos_cutted())
		{
			for (auto const process : set)
			{
				if (cut_required(process, cut_result))
				{
					borns.set(process, T());
				}
			}
		}

		auto const pdfs = luminosities_.pdfs(info.x1(), info.x2(),
			scales.factorization());

		T result = fold(pdfs, borns, set);
		result *= T(0.5) / info.energy_squared();
		result *= conversion_constant_;

		return result;
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
	L luminosities_;
	S scale_setter_;
	T conversion_constant_;

	T old_renormalization_scale_;
};

template <class T, class M, class C, class R, class L, class S>
using observables_born_like_type = observables_born_like<T,
	typename std::decay<M>::type, typename std::decay<C>::type,
	typename std::decay<R>::type, typename std::decay<L>::type,
	typename std::decay<S>::type>;

template <class T, class M, class C, class R, class L, class S>
inline observables_born_like_type <T, M, C, R, L, S> make_observables_born_like(
	M&& matrix_elements,
	C&& cuts,
	R&& recombiner,
	L&& luminosities,
	S&& scale_setter,
	T conversion_constant
) {
	return observables_born_like_type<T, M, C, R, L, S>(
		std::forward<M>(matrix_elements),
		std::forward<C>(cuts),
		std::forward<R>(recombiner),
		std::forward<L>(luminosities),
		std::forward<S>(scale_setter),
		conversion_constant
	);
}

}

#endif
