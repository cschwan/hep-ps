#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"
#include "hep/mc/projector.hpp"

#include "hep/ps/ab_terms.hpp"
#include "hep/ps/correction_type.hpp"
#include "hep/ps/cut_result.hpp"
#include "hep/ps/final_state.hpp"
#include "hep/ps/fini_integrand.hpp"
#include "hep/ps/finite_parts.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/parton_type.hpp"
#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scales.hpp"
#include "hep/ps/trivial_cutter.hpp"
#include "hep/ps/trivial_recombiner.hpp"

#include "test_phase_space_generator.hpp"

#include <catch.hpp>

#include <cstddef>
#include <functional>
#include <vector>

using T = HEP_TYPE_T;

std::vector<hep::scales<T>> global_scales;

template <typename T>
class test_fini_class
{
public:
	test_fini_class(
		std::size_t alphas_power,
		T alphas,
		bool dynamic_scale
	)
		: alphas_power_{alphas_power}
		, alphas_{alphas}
		, dynamic_scale_{dynamic_scale}
	{
	}

	// SUBTRACTION MEMBER FUNCTIONS

	void insertion_terms(
		hep::insertion_term const& term,
		std::vector<hep::scales<T>> const& scales,
		std::vector<T> const& /*phase_space*/,
		T /*x*/,
		T eta,
		std::vector<hep::ab_terms<T>>& results
	) const {
		results.clear();

		CHECK( term.type() == hep::insertion_term_type::born );

		for (auto const& scale : scales)
		{
			constexpr auto quark = hep::parton_type::quark;
			hep::ab_terms<T> term;
			term.a[quark][quark] = T(1.0);
			term.b[quark][quark] = (eta - T(1.0)) * T(2.0) + T(1.0);

			results.push_back(term);
		}
	}

	void insertion_terms2(
		hep::insertion_term const& /*term*/,
		std::vector<hep::scales<T>> const& /*scales*/,
		std::vector<T> const& /*phase_space*/,
		std::vector<T>& /*results*/
	) const {
		// TODO: test this call
	}

	// DISTRIBUTION MEMBER FUNCTIONS

	template <typename I>
	void operator()(
		std::vector<T> const& /*phase_space*/,
		T /*rapidity_shift*/,
		hep::cut_result_with_info<I> const&,
		std::vector<hep::neg_pos_results<T>> const& results,
		std::vector<hep::neg_pos_results<T>> const& pdf_results,
		hep::projector<T>&
	) {
		using std::pow;

		CHECK( results.size() == global_scales.size() );

		std::vector<T> alphas;
		eval_alphas(global_scales, alphas);

		for (std::size_t i = 0; i != results.size(); ++i)
		{
			T const muf = global_scales.at(i).factorization();
			T const mur = global_scales.at(i).renormalization();
			T const couplings = pow(alphas.at(i), T(alphas_power_));
			T const pdfa = muf;
			T const pdfb = muf;
			T const fini = T(2.0) * couplings;

			T const ref_result = pdfa * pdfb * fini;

			CHECK( results.at(i).neg == ref_result );
			CHECK( results.at(i).pos == ref_result );
		}

		for (std::size_t i = 0; i != pdf_results.size(); ++i)
		{
			T const muf = global_scales.front().factorization();
			T const mur = global_scales.front().renormalization();
			T const couplings = pow(alphas.front(), T(alphas_power_));
			T const pdfa = T(i+1) * muf;
			T const pdfb = T(i+1) * muf;
			T const fini = T(2.0) * couplings;

			T const ref_result = pdfa * pdfb * fini;

			CHECK( pdf_results.at(i).neg == ref_result );
			CHECK( pdf_results.at(i).pos == ref_result );
		}
	}

	// MATRIX ELEMENT MEMBER FUNCTIONS

	void correlated_me(
		std::vector<T> const& /*phase_space*/,
		hep::initial_state_set set,
		std::vector<hep::initial_state_map<T>>& results
	) {
		for (auto const state : set)
		{
			if (parton_type_of(state_parton_one(state)) !=
				parton_type_of(state_parton_two(state)))
			{
				continue;
			}

			results.at(0).emplace_back(state,
				pow(alphas_, T(alphas_power_)) * T(1.0));
			results.at(1).emplace_back(state,
				pow(alphas_, T(alphas_power_)) * T(1.0));
		}
	}

	void alphas(T alphas)
	{
		if (!dynamic_scale_)
		{
			CHECK( alphas_ == alphas );
		}
	}

	std::size_t alphas_power() const
	{
		return alphas_power_;
	}

	std::vector<hep::insertion_term> insertion_terms() const
	{
		return { hep::insertion_term{0, hep::correction_type::ew},
			hep::insertion_term{1, hep::correction_type::ew} };
	}

	std::vector<hep::final_state> final_states() const
	{
		return {
			hep::final_state::quark_gluon,
			hep::final_state::quark_gluon
		};
	}

	// SCALE SETTER MEMBER FUNCTIONS

	void operator()(
		std::vector<T> const&,
		std::vector<hep::scales<T>>& scales
	) {
		if (dynamic_scale_)
		{
			for (auto& scale : global_scales)
			{
				T const muf = T(2.0) * scale.factorization();
				T const mu = scale.regularization();
				T const mur = T(2.0) * scale.renormalization();

				scale = hep::scales<T>{muf, mu, mur};
			}
		}

		scales.assign(global_scales.begin(), global_scales.end());
	}

	bool dynamic() const
	{
		return dynamic_scale_;
	}

	// PDF MEMBER FUNCTIONS

	void eval_alphas(
		std::vector<hep::scales<T>> const& scales,
		std::vector<T>& alphas
	) {
		CHECK( scales.size() == global_scales.size() );

		for (std::size_t i = 0; i != scales.size(); ++i)
		{
			CHECK( scales.at(i).factorization() ==
				global_scales.at(i).factorization() );
			CHECK( scales.at(i).renormalization() ==
				global_scales.at(i).renormalization() );

			alphas.push_back(alphas_ * scales.at(i).renormalization() /
				global_scales.front().renormalization());
		}
	}

	std::size_t count() const
	{
		return 10;
	}

	void eval(
		T x,
		std::vector<hep::scales<T>> const& scales,
		std::vector<hep::parton_array<T>>& scale_pdfs,
		std::vector<hep::parton_array<T>>& uncertainty_pdfs
	) {
		CHECK( x >= T{} );
		CHECK( x <= T(1.0) );
		CHECK( scales.size() == global_scales.size() );

		scale_pdfs.clear();
		uncertainty_pdfs.clear();
		scale_pdfs.resize(scales.size());
		uncertainty_pdfs.resize(count());

		for (std::size_t i = 0; i != scales.size(); ++i)
		{
			for (auto const parton : hep::parton_list())
			{
				scale_pdfs.at(i)[parton] = scales.at(i).factorization();
			}
		}

		for (std::size_t i = 0; i != count(); ++i)
		{
			for (auto const parton : hep::parton_list())
			{
				uncertainty_pdfs.at(i)[parton] =
					scales.front().factorization() * T(i+1);
			}
		}
	}

	void register_partons(hep::parton_set)
	{
	}

private:
	std::size_t alphas_power_;
	T alphas_;
	bool dynamic_scale_;
};

void test_fini_integrand(
	hep::initial_state_set set,
	std::size_t alphas_power,
	T alphas,
	bool dynamic_scale
) {
	CAPTURE( alphas_power );
	CAPTURE( alphas );
	CAPTURE( dynamic_scale );

	// reset scales
	global_scales = {
		hep::scales<T>{         T(10.0), T(10.0),          T(10.0)},
		hep::scales<T>{T(0.5) * T(10.0), T(10.0),          T(10.0)},
		hep::scales<T>{         T(10.0), T(10.0), T(0.5) * T(10.0)},
		hep::scales<T>{T(2.0) * T(10.0), T(10.0),          T(10.0)},
		hep::scales<T>{         T(10.0), T(10.0), T(2.0) * T(10.0)},
		hep::scales<T>{T(0.5) * T(10.0), T(10.0), T(0.5) * T(10.0)},
		hep::scales<T>{T(2.0) * T(10.0), T(10.0), T(2.0) * T(10.0)}
	};

	// number of final states is not really interesting here
	test_phase_space_generator<T> generator{2, 1};
	test_fini_class<T> fini{alphas_power, alphas, dynamic_scale};

	// born gets copied, therefore `global_scales` must be global
	auto integrand = hep::make_fini_integrand<T>(
		fini,
		fini,
		hep::trivial_cutter<T>{},
		hep::trivial_recombiner<T>{},
		fini,
		fini,
		fini,
		set,
		T(1.0),
		hep::finite_parts::insertion_term
	);

	hep::multi_channel(
		hep::make_multi_channel_integrand<T>(
			std::ref(*integrand),
			generator.dimensions(),
			std::ref(generator),
			generator.map_dimensions(),
			generator.channels(),
			std::vector<hep::distribution_parameters<T>>{}
		),
		std::vector<std::size_t>{10}
	);
}

TEST_CASE("fini integrand static scale", "[fini_integrand]")
{
	// choose an initial state that doesn't have a symmetry factor
	hep::initial_state_set set{hep::initial_state::cq_uq};

	test_fini_integrand(set, 0, T(10.0), false);
	test_fini_integrand(set, 1, T(10.0), false);
	test_fini_integrand(set, 2, T(10.0), false);
	test_fini_integrand(set, 3, T(10.0), false);
	test_fini_integrand(set, 4, T(10.0), false);
}

TEST_CASE("fini integrand dynamic scale", "[fini_integrand]")
{
	// choose an initial state that doesn't have a symmetry factor
	hep::initial_state_set set{hep::initial_state::cq_uq};

	test_fini_integrand(set, 0, T(10.0), true);
	test_fini_integrand(set, 1, T(10.0), true);
	test_fini_integrand(set, 2, T(10.0), true);
	test_fini_integrand(set, 3, T(10.0), true);
	test_fini_integrand(set, 4, T(10.0), true);
}
