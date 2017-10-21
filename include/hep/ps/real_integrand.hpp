#ifndef HEP_PS_REAL_INTEGRAND_HPP
#define HEP_PS_REAL_INTEGRAND_HPP

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
#include "hep/ps/non_zero_dipole.hpp"
#include "hep/ps/particle_type.hpp"
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

namespace
{

inline void adjust_indices(
	std::vector<std::size_t> const& indices,
	std::size_t unresolved,
	std::vector<std::size_t>& result
) {
	result = indices;

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
}

}

namespace hep
{

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
class real_integrand : public ps_integrand<T>
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
	real_integrand(
		MatrixElements&& matrix_elements,
		Subtraction&& subtraction,
		Cuts&& cuts,
		Recombiner&& recombiner,
		Pdfs&& pdfs,
		ScaleSetter&& scale_setter,
		Distributions&& distributions,
		initial_state_set set,
		T hbarc2,
		T alpha_min
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
		, alpha_min_(alpha_min)
		, pdfsx1_(pdfs_.count())
		, pdfsx2_(pdfs_.count())
	{
		results_.reserve(pdfs_.count());
		dipole_recombination_candidates_.reserve(
			matrix_elements_.real_recombination_candidates().size());
		non_zero_dipoles_.reserve(matrix_elements_.dipoles().size());
	}

	T eval(
		std::vector<T> const& real_phase_space,
		luminosity_info<T> const& info,
		hep::projector<T>& projector
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

		auto set = set_;

		non_zero_dipoles_.clear();

		for (auto const dipole_with_set : matrix_elements_.dipoles())
		{
			auto const& dipole = dipole_with_set.dipole();

			std::vector<T> phase_space(real_phase_space.size() - 4);

			// map the real phase space on the dipole phase space
			auto const invariants = subtraction_.map_phase_space(
				real_phase_space, phase_space, dipole);

			if (invariants.alpha < alpha_min_)
			{
				// remove all initial states from the set that share this dipole
				set.subtract(dipole_with_set.set());
				continue;
			}

			adjust_indices(
				matrix_elements_.real_recombination_candidates(),
				dipole.unresolved(),
				dipole_recombination_candidates_
			);

			auto const dipole_recombined = recombiner_.recombine(
				phase_space,
				phase_space,
				dipole_recombination_candidates_,
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

			non_zero_dipoles_.emplace_back(
				std::move(phase_space),
				invariants,
				dipole,
				dipole_cut_result
			);
		}

		// if there are neither dipoles nor real matrix elements stop here
		if (set.empty() || (non_zero_dipoles_.empty() &&
			real_cut_result.neg_cutted() && real_cut_result.pos_cutted()))
		{
			return T();
		}

		auto const scales = scale_setter_(recombined_real_phase_space);

		// only set renormalization scale if it changed
		if (scales.renormalization() != old_renormalization_scale_)
		{
			old_renormalization_scale_ = scales.renormalization();

			T const alphas = pdfs_.alphas(scales.renormalization());
			matrix_elements_.parameters(scales.renormalization(), alphas);
		}

		pdfs_.eval(info.x1(), scales.factorization(), pdfsx1_);
		pdfs_.eval(info.x2(), scales.factorization(), pdfsx2_);

		T const factor = T(0.5) * hbarc2_ / info.energy_squared();

		neg_pos_results<T> result;

		for (auto const non_zero_dipole : non_zero_dipoles_)
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

			std::size_t size = pdfsx1_.size();
			results_.clear();

			for (std::size_t pdf = 0; pdf != size; ++pdf)
			{
				results_.push_back(convolute(
					pdfsx1_.at(pdf),
					pdfsx2_.at(pdf),
					dipole_me,
					set,
					-function * factor,
					dipole_cut_result
				));
			}

			distributions_(phase_space, dipole_cut_result, results_, shift,
				event_type::born_like_n, projector);

			result += results_.front();
		}

		if (!real_cut_result.neg_cutted() || !real_cut_result.pos_cutted())
		{
			auto const reals = matrix_elements_.reals(real_phase_space, set);

			std::size_t const size = pdfsx1_.size();
			results_.clear();

			for (std::size_t pdf = 0; pdf != size; ++pdf)
			{
				results_.push_back(convolute(
					pdfsx1_.at(pdf),
					pdfsx2_.at(pdf),
					reals,
					set,
					factor,
					real_cut_result
				));
			}

			result += results_.front();

			distributions_(recombined_real_phase_space, real_cut_result,
				results_, shift, event, projector);
		}

		return result.neg + result.pos;
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
	T alpha_min_;

	T old_renormalization_scale_;

	using info_type = typename decltype (cuts_.cut(std::vector<T>(), T{},
		event_type{}))::info_t;

	std::vector<parton_array<T>> pdfsx1_;
	std::vector<parton_array<T>> pdfsx2_;
	std::vector<neg_pos_results<T>> results_;
	std::vector<std::size_t> dipole_recombination_candidates_;
	std::vector<non_zero_dipole<T, info_type>> non_zero_dipoles_;
};

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
using real_integrand_type = real_integrand<T, std::decay_t<M>, std::decay_t<S>,
	std::decay_t<C>, std::decay_t<R>, std::decay_t<P>, std::decay_t<U>,
	std::decay_t<D>>;

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
inline std::unique_ptr<ps_integrand<T>> make_real_integrand(
	M&& matrix_elements,
	S&& subtraction,
	C&& cuts,
	R&& recombiner,
	P&& pdfs,
	U&& scale_setter,
	D&& distributions,
	initial_state_set set,
	T hbarc2,
	T alpha_min = T()
) {
	return std::make_unique<real_integrand_type<T, M, S, C, R, P, U, D>>(
		std::forward<M>(matrix_elements),
		std::forward<S>(subtraction),
		std::forward<C>(cuts),
		std::forward<R>(recombiner),
		std::forward<P>(pdfs),
		std::forward<U>(scale_setter),
		std::forward<D>(distributions),
		set,
		hbarc2,
		alpha_min
	);
}

}

#endif
