#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"
#include "hep/mc/projector.hpp"

#include "hep/ps/correction_type.hpp"
#include "hep/ps/cut_result.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_invariants.hpp"
#include "hep/ps/dipole_type.hpp"
#include "hep/ps/final_state.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/neg_pos_results.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/real_integrand.hpp"
#include "hep/ps/recombined_state.hpp"
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
class test_real_class
{
public:
	test_real_class(
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

	hep::dipole_invariants<T> map_phase_space(
		std::vector<T> const& /*real_phase_space*/,
		std::vector<T>& /*born_phase_space*/,
		hep::dipole const& /*dipole_info*/
	) {
		return {};
	}

	T fermion_function(
		hep::dipole const& /*dipole_info*/,
		hep::dipole_invariants<T> const& /*invariants*/
	) {
		return T(1.0);
	}

	// DISTRIBUTION MEMBER FUNCTIONS

	template <typename I>
	void operator()(
		std::vector<T> const& phase_space,
		T /*rapidity_shift*/,
		hep::cut_result_with_info<I> const& /*cut_result*/,
		std::vector<hep::neg_pos_results<T>> const& results,
		std::vector<hep::neg_pos_results<T>> const& pdf_results,
		hep::projector<T>&
	) {
		using std::pow;

		CHECK( results.size() == global_scales.size() );

		std::vector<T> alphas;
		eval_alphas(global_scales, alphas);

		T factor;

		switch (phase_space.size())
		{
		case 20:
			// real matrix element
			factor = T(2.0);
			break;

		case 16:
			// dipole
			factor = T(-1.0);
			break;

		default:
			assert( false );
		}

		for (std::size_t i = 0; i != results.size(); ++i)
		{
			T const muf = global_scales.at(i).factorization();
			T const mur = global_scales.at(i).renormalization();
			T const couplings = pow(alphas.at(i), T(alphas_power_));
			T const pdfa = muf;
			T const pdfb = muf;
			T const real = factor * couplings;

			T const ref_result = pdfa * pdfb * real;

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
			T const real = factor * couplings;

			T const ref_result = pdfa * pdfb * real;

			CHECK( pdf_results.at(i).neg == ref_result );
			CHECK( pdf_results.at(i).pos == ref_result );
		}
	}

	// MATRIX ELEMENT MEMBER FUNCTIONS

	void reals(
		std::vector<T> const&,
		hep::initial_state_set set,
		std::vector<hep::scales<T>> const& scales,
		std::vector<hep::initial_state_map<T>>& results
	) const {
		using std::pow;

		hep::initial_state_map<T> result;

		for (auto const state : set)
		{
			result.emplace_back(state, pow(alphas_, T(alphas_power_)) * T(2.0));
		}

		results.assign(scales.size(), result);
	}

	void dipole_me(
		hep::dipole const& /*dipole*/,
		std::vector<T> const& /*phase_space*/,
		hep::initial_state_set set,
		std::vector<hep::scales<T>> const& scales,
		std::vector<hep::initial_state_map<T>>& results
	) {
		hep::initial_state_map<T> result;

		for (auto const state : set)
		{
			result.emplace_back(state, pow(alphas_, T(alphas_power_)) * T(1.0));
		}

		results.assign(scales.size(), result);
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

	std::vector<hep::dipole> dipoles() const
	{
		return { hep::dipole{ 2, 4, 3, hep::particle_type::fermion,
			hep::particle_type::boson, hep::particle_type::fermion,
			hep::correction_type::qcd} };
	}

	std::vector<hep::final_state> final_states() const
	{
		return {
			hep::final_state::quark_gluon,
			hep::final_state::quark_gluon
		};
	}

	std::vector<hep::final_state> final_states_real() const
	{
		return {
			hep::final_state::quark_gluon,
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
		CHECK( x < T(1.0) );
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

void test_real_integrand(
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
	test_phase_space_generator<T> generator{2 + 1};
	test_real_class<T> real{alphas_power, alphas, dynamic_scale};

	// born gets copied, therefore `global_scales` must be global
	auto integrand = hep::make_real_integrand<T>(
		real,
		real,
		hep::trivial_cutter<T>{},
		hep::trivial_recombiner<T>{},
		real,
		real,
		real,
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

TEST_CASE("real integrand static scale", "[real_integrand]")
{
	// choose an initial state that doesn't have a symmetry factor
	hep::initial_state_set set{hep::initial_state::cq_uq};

	test_real_integrand(set, 0, T(10.0), false);
	test_real_integrand(set, 1, T(10.0), false);
	test_real_integrand(set, 2, T(10.0), false);
	test_real_integrand(set, 3, T(10.0), false);
	test_real_integrand(set, 4, T(10.0), false);
}

TEST_CASE("real integrand dynamic scale", "[real_integrand]")
{
	// choose an initial state that doesn't have a symmetry factor
	hep::initial_state_set set{hep::initial_state::cq_uq};

	test_real_integrand(set, 0, T(10.0), true);
	test_real_integrand(set, 1, T(10.0), true);
	test_real_integrand(set, 2, T(10.0), true);
	test_real_integrand(set, 3, T(10.0), true);
	test_real_integrand(set, 4, T(10.0), true);
}
