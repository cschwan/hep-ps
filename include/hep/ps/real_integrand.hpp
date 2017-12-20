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
#include "hep/ps/final_state.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/luminosity_info.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/non_zero_dipole.hpp"
#include "hep/ps/particle_type.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/ps_integrand.hpp"
#include "hep/ps/recombined_state.hpp"
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
		, final_states_real_(matrix_elements_.final_states_real())
		, final_states_dipole_(matrix_elements_.final_states())
		, pdf_pdfsx1_(pdfs_.count())
		, pdf_pdfsx2_(pdfs_.count())
		, alphas_power_(matrix_elements_.alphas_power())
	{
		std::size_t const fs = final_states_real_.size();

		recombined_ps_.reserve(4 * (fs + 2));
		recombined_states_.reserve(fs);
		recombined_dipole_states_.reserve(fs - 1);
		pdf_results_.reserve(pdfs_.count());

		non_zero_dipoles_.reserve(matrix_elements_.dipoles().size());

		if (!scale_setter_.dynamic())
		{
			set_scales(std::vector<T>());
		}
	}

	T eval(
		std::vector<T> const& real_phase_space,
		luminosity_info<T> const& info,
		hep::projector<T>& projector
	) override {
		recombiner_.recombine(
			real_phase_space,
			final_states_real_,
			recombined_ps_,
			recombined_states_
		);

		auto const real_cut_result = cuts_.cut(
			recombined_ps_,
			info.rapidity_shift(),
			recombined_states_
		);

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

			recombiner_.recombine(
				phase_space,
				final_states_dipole_,
				phase_space,
				recombined_dipole_states_
			);

			auto const dipole_cut_result = cuts_.cut(
				phase_space,
				info.rapidity_shift(),
				recombined_dipole_states_
			);

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

		if (scale_setter_.dynamic())
		{
			set_scales(recombined_ps_);
		}

		pdfsx1_.clear();
		pdfsx1_.resize(scales_.size());
		pdfsx2_.clear();
		pdfsx2_.resize(scales_.size());

		pdfs_.eval(info.x1(), scales_, pdfsx1_);
		pdfs_.eval(info.x2(), scales_, pdfsx2_);

		assert( pdfsx1_.size() == scales_.size() );
		assert( pdfsx2_.size() == scales_.size() );

		if (pdfs_.count() > 1)
		{
			pdf_pdfsx1_.clear();
			pdf_pdfsx1_.resize(pdfs_.count());
			pdf_pdfsx2_.clear();
			pdf_pdfsx2_.resize(pdfs_.count());

			T const muf = scales_.front().factorization();

			pdfs_.eval(info.x1(), muf, pdf_pdfsx1_);
			pdfs_.eval(info.x2(), muf, pdf_pdfsx2_);

			assert( pdf_pdfsx1_.size() == pdfs_.count() );
			assert( pdf_pdfsx2_.size() == pdfs_.count() );
		}

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

			results_.clear();

			for (std::size_t i = 0; i != scales_.size(); ++i)
			{
				results_.push_back(convolute(
					pdfsx1_.at(i),
					pdfsx2_.at(i),
					dipole_me,
					set,
					-function * factors_.at(i) * factor,
					dipole_cut_result
				));
			}

			if (pdfs_.count() > 1)
			{
				pdf_results_.clear();

				for (std::size_t pdf = 0; pdf != pdfs_.count(); ++pdf)
				{
					pdf_results_.push_back(convolute(
						pdf_pdfsx1_.at(pdf),
						pdf_pdfsx2_.at(pdf),
						dipole_me,
						set,
						-function * factor,
						dipole_cut_result
					));
				}
			}

			distributions_(
				phase_space,
				info.rapidity_shift(),
				dipole_cut_result,
				results_,
				pdf_results_,
				projector
			);

			result += results_.front();
		}

		if (!real_cut_result.neg_cutted() || !real_cut_result.pos_cutted())
		{
			auto const reals = matrix_elements_.reals(real_phase_space, set);

			results_.clear();

			for (std::size_t i = 0; i != scales_.size(); ++i)
			{
				results_.push_back(convolute(
					pdfsx1_.at(i),
					pdfsx2_.at(i),
					reals,
					set,
					factors_.at(i) * factor,
					real_cut_result
				));
			}

			if (pdfs_.count() > 1)
			{
				pdf_results_.clear();

				for (std::size_t pdf = 0; pdf != pdfs_.count(); ++pdf)
				{
					pdf_results_.push_back(convolute(
						pdf_pdfsx1_.at(pdf),
						pdf_pdfsx2_.at(pdf),
						reals,
						set,
						factor,
						real_cut_result
					));
				}
			}

			distributions_(
				recombined_ps_,
				info.rapidity_shift(),
				real_cut_result,
				results_,
				pdf_results_,
				projector
			);

			result += results_.front();
		}

		return result.neg + result.pos;
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
	S subtraction_;
	C cuts_;
	R recombiner_;
	P pdfs_;
	U scale_setter_;
	D distributions_;
	initial_state_set set_;
	T hbarc2_;
	T alpha_min_;

	using info_type = typename decltype (cuts_.cut(std::vector<T>(), T{},
		std::vector<recombined_state>{}))::info_t;

	std::vector<T> recombined_ps_;
	std::vector<recombined_state> recombined_states_;
	std::vector<recombined_state> recombined_dipole_states_;
	std::vector<final_state> final_states_real_;
	std::vector<final_state> final_states_dipole_;
	std::vector<parton_array<T>> pdfsx1_;
	std::vector<parton_array<T>> pdfsx2_;
	std::vector<parton_array<T>> pdf_pdfsx1_;
	std::vector<parton_array<T>> pdf_pdfsx2_;
	std::vector<neg_pos_results<T>> pdf_results_;
	std::vector<neg_pos_results<T>> results_;
	std::vector<scales<T>> scales_;
	std::vector<T> factors_;
	std::vector<non_zero_dipole<T, info_type>> non_zero_dipoles_;
	T alphas_power_;
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
