#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"
#include "hep/mc/projector.hpp"

#include "hep/ps/abc_terms.hpp"
#include "hep/ps/cut_result.hpp"
#include "hep/ps/event_type.hpp"
#include "hep/ps/fini_integrand.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/scales.hpp"
#include "hep/ps/trivial_cutter.hpp"
#include "hep/ps/trivial_recombiner.hpp"

#include "catch.hpp"
#include "test_phase_space_generator.hpp"

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
		T /*eta*/,
		std::vector<hep::abc_terms<T>>& results
	) const {
		CHECK( term.type() == hep::insertion_term_type::born );

		for (auto const& scale : scales)
		{
			hep::abc_terms<T> term;
			term.a[hep::parton_type::quark][hep::parton_type::quark] = T(1.0);
			term.b[hep::parton_type::quark][hep::parton_type::quark] = T(2.0);
			term.c[hep::parton_type::quark][hep::parton_type::quark] = T(1.0);

			results.push_back(term);
		}
	}

	void insertion_terms2(
		hep::insertion_term const& term,
		std::vector<hep::scales<T>> const& scales,
		std::vector<T> const& /*phase_space*/,
		std::vector<T>& results
	) const {
		CHECK( term.type() == hep::insertion_term_type::born );

		for (auto const& scale : scales)
		{
			// TODO: test this call
			results.emplace_back();
		}
	}

	// DISTRIBUTION MEMBER FUNCTIONS

	template <typename I>
	void operator()(
		std::vector<T> const& /*phase_space*/,
		hep::cut_result_with_info<I> const&,
		std::vector<hep::neg_pos_results<T>> const& results,
		std::vector<hep::neg_pos_results<T>> const& pdf_results,
		T,
		hep::event_type event_type,
		hep::projector<T>&
	) {
		using std::pow;

		CHECK( results.size() == global_scales.size() );
		CHECK( event_type == hep::event_type::born_like_n );

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

	std::array<hep::initial_state_array<T>, 1> correlated_me(
		std::vector<T> const& /*phase_space*/,
		hep::initial_state_set set
	) {
		hep::initial_state_array<T> result;

		for (auto const state : set)
		{
			result[state] = pow(alphas_, T(alphas_power_)) * T(1.0);
		}

		return { result };
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

	std::array<hep::insertion_term, 1> insertion_terms() const
	{
		return { hep::insertion_term{hep::insertion_term_type::born} };
	}

	std::vector<std::size_t> born_recombination_candidates() const
	{
		return { 2, 3 };
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
				T const mur = T(2.0) * scale.renormalization();

				scale = hep::scales<T>{muf, mur};
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
		std::vector<hep::parton_array<T>>& pdfs
	) {
		CHECK( x >= T{} );
		CHECK( x <= T(1.0) );
		CHECK( scales.size() == global_scales.size() );
		CHECK( pdfs.size() == scales.size() );

		for (std::size_t i = 0; i != scales.size(); ++i)
		{
			for (auto const parton : hep::parton_list())
			{
				pdfs.at(i)[parton] = scales.at(i).factorization();
			}
		}
	}

	void eval(T x, T scale, std::vector<hep::parton_array<T>>& pdfs)
	{
		CHECK( x >= T{} );
		CHECK( x <= T(1.0) );
		CHECK( scale > T{} );
		CHECK( pdfs.size() == count() );

		for (std::size_t i = 0; i != count(); ++i)
		{
			for (auto const parton : hep::parton_list())
			{
				pdfs.at(i)[parton] = scale * T(i+1);
			}
		}
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
		hep::scales<T>{         T(10.0),          T(10.0)},
		hep::scales<T>{T(0.5) * T(10.0),          T(10.0)},
		hep::scales<T>{         T(10.0), T(0.5) * T(10.0)},
		hep::scales<T>{T(2.0) * T(10.0),          T(10.0)},
		hep::scales<T>{         T(10.0), T(2.0) * T(10.0)},
		hep::scales<T>{T(0.5) * T(10.0), T(0.5) * T(10.0)},
		hep::scales<T>{T(2.0) * T(10.0), T(2.0) * T(10.0)}
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
		T(1.0)
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
	hep::initial_state_set set{hep::initial_state::q43_cu};

	test_fini_integrand(set, 0, T(10.0), false);
	test_fini_integrand(set, 1, T(10.0), false);
	test_fini_integrand(set, 2, T(10.0), false);
	test_fini_integrand(set, 3, T(10.0), false);
	test_fini_integrand(set, 4, T(10.0), false);
}

TEST_CASE("fini integrand dynamic scale", "[fini_integrand]")
{
	// choose an initial state that doesn't have a symmetry factor
	hep::initial_state_set set{hep::initial_state::q43_cu};

	test_fini_integrand(set, 0, T(10.0), true);
	test_fini_integrand(set, 1, T(10.0), true);
	test_fini_integrand(set, 2, T(10.0), true);
	test_fini_integrand(set, 3, T(10.0), true);
	test_fini_integrand(set, 4, T(10.0), true);
}
