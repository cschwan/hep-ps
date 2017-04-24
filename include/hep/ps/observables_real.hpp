#ifndef HEP_PS_OBSERVABLES_REAL_HPP
#define HEP_PS_OBSERVABLES_REAL_HPP

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

#include "hep/ps/distributions.hpp"
#include "hep/ps/event_type.hpp"
#include "hep/ps/fold.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/non_zero_dipole.hpp"
#include "hep/ps/observables.hpp"
#include "hep/ps/particle_type.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

inline std::vector<std::size_t> adjust_indices(
	std::vector<std::size_t> const& indices,
	std::size_t unresolved
) {
	std::vector<std::size_t> result(indices);

	auto const begin = std::find(result.begin(), result.end(), unresolved);
	auto end = result.end();

	assert( begin != end );

	auto next = std::next(begin);

	// decrease all indices following the index for the unresolved by one
	std::transform(next, end, next, [](std::size_t v) { return v - 1; });
	// rotate the unresolved index to the end of the vector
	std::rotate(begin, next, end);
	// remove the unresolved index
	result.erase(std::prev(end));

	return result;
}

}

namespace hep
{

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
class observables_real : public observables<T>
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
	observables_real(
		MatrixElements&& matrix_elements,
		Subtraction&& subtraction,
		Cuts&& cuts,
		Recombiner&& recombiner,
		Pdf&& pdf,
		ScaleSetter&& scale_setter,
		Distributions&& distributions,
		T hbarc2,
		T alpha_min
	)
		: matrix_elements_(std::forward<MatrixElements>(matrix_elements))
		, subtraction_(std::forward<Subtraction>(subtraction))
		, cuts_(std::forward<Cuts>(cuts))
		, recombiner_(std::forward<Recombiner>(recombiner))
		, pdf_(std::forward<Pdf>(pdf))
		, scale_setter_(std::forward<ScaleSetter>(scale_setter))
		, distributions_(std::forward<Distributions>(distributions))
		, hbarc2_(hbarc2)
		, alpha_min_(alpha_min)
	{
		static_assert (std::is_base_of<hep::distributions<T>, D>::value,
			"`D` must be a type deriving from `hep::distributions<T>`");
	}

	T eval(
		std::vector<T> const& real_phase_space,
		luminosity_info<T> const& info,
		initial_state_set set
	) override {
		std::vector<T> recombined_real_phase_space(real_phase_space.size());

		auto const recombined = recombiner_.recombine(
			real_phase_space,
			recombined_real_phase_space,
			matrix_elements_.real_recombination_candidates(),
			1
		);

		event_type event;

		switch (recombined)
		{
		case 0: event = event_type::inclusive_n_plus_1; break;
		case 1: event = event_type::born_like_n; break;
		default:
			event = event_type::other;
		}

		T const shift = info.rapidity_shift();

		using cut_result_type = decltype (cuts_.cut(recombined_real_phase_space,
			shift, event));

		cut_result_type real_cut_result;

		if (event == event_type::inclusive_n_plus_1 ||
			event == event_type::born_like_n)
		{
			real_cut_result = cuts_.cut(recombined_real_phase_space, shift,
				event);
		}

		using info_type = typename cut_result_type::info_t;

		std::vector<non_zero_dipole<T, info_type>> non_zero_dipoles;

		for (auto const dipole_with_set : matrix_elements_.dipoles())
		{
			auto const& dipole = dipole_with_set.dipole();

			std::vector<T> phase_space(real_phase_space.size() - 4);

			// map the real phase space on the dipole phase space
			auto const invariants = subtraction_.map_phase_space(
				real_phase_space, phase_space, dipole);

			if (invariants.adipole < alpha_min_)
			{
				// remove all initial states from the set that share this dipole
				set.subtract(dipole_with_set.set());
				continue;
			}

			auto const dipole_recombination_candidates = adjust_indices(
				matrix_elements_.real_recombination_candidates(),
				dipole.unresolved());

			auto const dipole_recombined = recombiner_.recombine(
				phase_space,
				phase_space,
				dipole_recombination_candidates,
				0
			);

			// check if it passed the recombination
			if (dipole_recombined > 0)
			{
				continue;
			}

			auto const dipole_cut_result = cuts_.cut(phase_space, shift,
				event_type::born_like_n);

			if (dipole_cut_result.neg_cutted() &&
				dipole_cut_result.pos_cutted())
			{
				continue;
			}

			non_zero_dipoles.emplace_back(
				std::move(phase_space),
				invariants,
				dipole,
				dipole_cut_result
			);
		}

		// if there are neither dipoles nor real matrix elements stop here
		if (set.empty() || (non_zero_dipoles.empty() &&
			real_cut_result.neg_cutted() && real_cut_result.pos_cutted()))
		{
			return T();
		}

		auto const scales = scale_setter_(recombined_real_phase_space);

		// only set renormalization scale if it changed
		if (scales.renormalization() != old_renormalization_scale_)
		{
			matrix_elements_.scale(scales.renormalization(), pdf_);
			old_renormalization_scale_ = scales.renormalization();
		}

		initial_state_array<T> reals;

		if (!real_cut_result.neg_cutted() || !real_cut_result.pos_cutted())
		{
			reals = matrix_elements_.reals(real_phase_space, set);
		}

		auto const pdfx1 = pdf_.pdf(info.x1(), scales.factorization());
		auto const pdfx2 = pdf_.pdf(info.x2(), scales.factorization());

		T const factor = T(0.5) * hbarc2_ / info.energy_squared();

		neg_pos_results<T> result;

		for (auto const non_zero_dipole : non_zero_dipoles)
		{
			auto const& dipole = non_zero_dipole.dipole();
			auto const& dipole_cut_result = non_zero_dipole.cut_result();
			auto const& phase_space = non_zero_dipole.phase_space();
			auto const& invariants = non_zero_dipole.invariants();

			bool const fermion_i = dipole.emitter_type() ==
				particle_type::fermion;
			bool const fermion_j = dipole.unresolved_type() ==
				particle_type::fermion;

			T function;

			if (fermion_i != fermion_j)
			{
				function = subtraction_.fermion_function(dipole,
					invariants);
			}
			else if (fermion_i && fermion_j)
			{
				// TODO: NYI
				assert( false );
			}
			else
			{
				// TODO: NYI
				assert( false );
			}

			auto const dipole_me = matrix_elements_.dipole_me(dipole,
				phase_space, set);
			auto const dipole_result = fold(pdfx1, pdfx2, dipole_me, set,
				-function * factor, dipole_cut_result);

			distributions_(phase_space, dipole_cut_result, dipole_result,
				shift, event_type::born_like_n);

			result += dipole_result;
		}

		auto const real_result = fold(pdfx1, pdfx2, reals, set, factor,
			real_cut_result);
		result += real_result;

		distributions_(recombined_real_phase_space, real_cut_result, real_result,
			shift, event);

		return result.neg + result.pos;
	}

	hep::distributions<T>& distributions() override
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
	T alpha_min_;

	T old_renormalization_scale_;
};

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
using observables_real_type = observables_real<T,
	typename std::decay<M>::type, typename std::decay<S>::type,
	typename std::decay<C>::type, typename std::decay<R>::type,
	typename std::decay<P>::type, typename std::decay<U>::type,
	typename std::decay<D>::type>;

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
inline std::unique_ptr<observables<T>> make_observables_real(
	M&& matrix_elements,
	S&& subtraction,
	C&& cuts,
	R&& recombiner,
	P&& pdf,
	U&& scale_setter,
	D&& distributions,
	T hbarc2,
	T alpha_min = T()
) {
	return std::unique_ptr<observables_real_type<T, M, S, C, R, P, U, D>>(
		new observables_real_type<T, M, S, C, R, P, U, D>(
			std::forward<M>(matrix_elements),
			std::forward<S>(subtraction),
			std::forward<C>(cuts),
			std::forward<R>(recombiner),
			std::forward<P>(pdf),
			std::forward<U>(scale_setter),
			std::forward<D>(distributions),
			hbarc2,
			alpha_min
	));
}

}

#endif
