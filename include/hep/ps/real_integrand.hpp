#ifndef HEP_PS_REAL_INTEGRAND_HPP
#define HEP_PS_REAL_INTEGRAND_HPP

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
#include "hep/ps/non_zero_dipole.hpp"
#include "hep/ps/particle_type.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/ps_integrand.hpp"
#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scales.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <memory>
#include <type_traits>
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
		, dipoles_(matrix_elements_.dipoles())
		, alphas_power_(matrix_elements_.alphas_power())
	{
		std::size_t const fs = final_states_real_.size();

		recombined_ps_.reserve(4 * (fs + 2));
		recombined_states_.reserve(fs);
		recombined_dipole_states_.reserve(fs - 1);
		pdf_results_.reserve((pdfs_.count() == 1) ? 0 : pdfs_.count());

		non_zero_dipoles_.reserve(dipoles_.size());
		dipole_phase_spaces_.resize(dipoles_.size());

		for (auto& phase_space : dipole_phase_spaces_)
		{
			phase_space.resize(4 * (fs + 1));
		}

		if (!scale_setter_.dynamic())
		{
			set_scales(std::vector<T>());
			results_.reserve(scales_.size());
			me_.resize(scales_.size());
			me_tmp_.resize(scales_.size());
		}

		pdfs_.register_partons(partons_in_initial_state_set(set));
	}

	T eval(
		std::vector<T> const& real_phase_space,
		luminosity_info<T> const& info,
		hep::projector<T>& projector
	) override {
		non_zero_dipoles_.clear();

		std::size_t index = 0;

		for (auto const& dipole : dipoles_)
		{
			auto& phase_space = dipole_phase_spaces_.at(index);

			// map the real phase space on the dipole phase space
			auto const invariants = subtraction_.map_phase_space(
				real_phase_space, phase_space, dipole);

			if (invariants.alpha < alpha_min_)
			{
				// if we apply a technical cut, completely throw away the point
				return T();
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
				index++,
				invariants,
				dipole,
				dipole_cut_result
			);
		}

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

		// if there are neither dipoles nor real matrix elements stop here
		if (non_zero_dipoles_.empty() && real_cut_result.neg_cutted() &&
			real_cut_result.pos_cutted())
		{
			return T();
		}

		if (scale_setter_.dynamic())
		{
			set_scales(recombined_ps_);
			results_.reserve(scales_.size());
			me_.resize(scales_.size());
			me_tmp_.resize(scales_.size());
		}

		pdfs_.eval(info.x1(), scales_, pdfsx1_, pdf_pdfsx1_);
		pdfs_.eval(info.x2(), scales_, pdfsx2_, pdf_pdfsx2_);

		assert( pdfsx1_.size() == scales_.size() );
		assert( pdfsx2_.size() == scales_.size() );
		assert( (pdfs_.count() == 1) || (pdf_pdfsx1_.size() == pdfs_.count()) );
		assert( (pdfs_.count() == 1) || (pdf_pdfsx2_.size() == pdfs_.count()) );

		T const factor = T(0.5) * hbarc2_ / info.energy_squared();

		neg_pos_results<T> result;

		for (auto const non_zero_dipole : non_zero_dipoles_)
		{
			auto const& dipole = non_zero_dipole.dipole();
			auto const& dipole_cut_result = non_zero_dipole.cut_result();
			auto const& invariants = non_zero_dipole.invariants();
			auto const& phase_space =
				dipole_phase_spaces_.at(non_zero_dipole.index());

			T function = T(1.0);

			if ((dipole.emitter_type() == particle_type::fermion) !=
				(dipole.unresolved_type() == particle_type::fermion))
			{
				function = -subtraction_.fermion_function(dipole, invariants);

				for (auto& me : me_)
				{
					me.clear();
				}

				matrix_elements_.dipole_me(
					dipole,
					phase_space,
					set_,
					scales_,
					me_
				);
			}
			else
			{
				std::array<std::array<T, 4>, 4> vectors;
				std::array<T, 4> functions;

				subtraction_.boson_function(dipole, invariants, phase_space,
					vectors, functions);

				for (std::size_t i = 0; i != vectors.size(); ++i)
				{
					for (auto& me : me_)
					{
						me.clear();
					}

					matrix_elements_.dipole_sc(
						dipole,
						phase_space,
						vectors[i],
						set_,
						scales_,
						(i == 0) ? me_ : me_tmp_
					);

					if (i == 0)
					{
						for (std::size_t j = 0; j != me_tmp_.size(); ++j)
						{
							for (std::size_t k = 0; me_tmp_[j].size(); ++k)
							{
								me_[j].emplace_back(me_tmp_[j][k].first,
									functions[i] * me_tmp_[j][k].second);
							}
						}
					}
					else
					{
						for (std::size_t j = 0; j != me_.size(); ++j)
						{
							for (std::size_t k = 0; me_[j].size(); ++k)
							{
								// order of `me_` and `me_tmp_` must be the same
								assert(me_[j][k].first ==
									me_tmp_[j].at(k).first);

								me_[j][k].second += functions[i] *
									me_tmp_.at(j).at(k).second;
							}
						}
					}
				}
			}

			convolute_mes_with_pdfs(
				results_,
				pdf_results_,
				pdfsx1_,
				pdfsx2_,
				pdf_pdfsx1_,
				pdf_pdfsx2_,
				me_,
				set_,
				factors_,
				function * factor,
				dipole_cut_result
			);

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
			for (auto& me : me_)
			{
				me.clear();
			}

			matrix_elements_.reals(real_phase_space, set_, scales_, me_);

			convolute_mes_with_pdfs(
				results_,
				pdf_results_,
				pdfsx1_,
				pdfsx2_,
				pdf_pdfsx1_,
				pdf_pdfsx2_,
				me_,
				set_,
				factors_,
				factor,
				real_cut_result
			);

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
	std::vector<std::vector<T>> dipole_phase_spaces_;
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
	std::vector<initial_state_map<T>> me_;
	std::vector<initial_state_map<T>> me_tmp_;
	std::vector<non_zero_dipole<T, info_type>> non_zero_dipoles_;
	std::vector<dipole> dipoles_;
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
