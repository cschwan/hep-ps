#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "hep/ps/mc_distributions.hpp"
#include "hep/ps/mc_observables.hpp"
#include "hep/ps/mc_phase_space_adapter.hpp"
#include "hep/ps/observables_born.hpp"
#include "hep/ps/observables_fini.hpp"
#include "hep/ps/observables_real.hpp"
#include "hep/ps/trivial_distributions.hpp"

#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"

#include "catch.hpp"
#include "test_structures.hpp"

#include <cstddef>
#include <functional>

using T = double;

hep::lusifer_constants<T> constants(
	T(125.09), T(4.0e-3),
	T(174.2), T(1.41),
	T(80.385), T(2.085),
	T(91.1876), T(2.4952)
);

TEST_CASE("test make_mc_observables_born", "[mc_observables]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};
	std::size_t count = 2;
	T const min_energy = T(10.0);
	T const cmf_energy = T(1000.0);

	auto generator = hep::mc_phase_space_adapter<T>(
		hep::make_lusifer_phase_space_generator<T>(
			min_energy,
			cmf_energy,
			"el~el mu mu~",
			constants
		)
	);

	hep::mc_observables<T> observables{hep::make_observables_born(
		test_matrix_elements<T>(set, count, 0, 0, 0),
		test_cuts<T>(count, false),
		test_recombiner<T>(count),
		test_pdf<T>(),
		test_scale_setter<T>(),
		hep::mc_trivial_distributions<T>(),
		set,
		T(1.0)
	)};

	auto const results = hep::multi_channel(
		hep::make_multi_channel_integrand<T>(
			std::ref(observables),
			generator.dimensions(),
			std::ref(generator),
			generator.map_dimensions(),
			generator.channels(),
			std::vector<hep::distribution_parameters<T>>()
		),
		std::vector<std::size_t>(1, 10)
	);

	CHECK( results.front().value() != T() );
	CHECK( results.front().error() != T() );
}

TEST_CASE("test make_mc_observables_fini", "[mc_observables]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};
	std::size_t count = 2;
	T const min_energy = T(10.0);
	T const cmf_energy = T(1000.0);

	auto generator = hep::mc_phase_space_adapter<T>(
		hep::make_lusifer_phase_space_generator<T>(
			min_energy,
			cmf_energy,
			"el~el mu mu~",
			constants,
			1
		)
	);

	hep::mc_observables<T> observables{hep::make_observables_fini(
		test_matrix_elements<T>(set, count, 0, 0, 0),
		test_subtraction<T>(0, 0, 0),
		test_cuts<T>(count, false),
		test_recombiner<T>(count),
		test_pdf<T>(),
		test_scale_setter<T>(),
		hep::mc_trivial_distributions<T>(),
		set,
		T(1.0),
		false
	)};

	auto const results = hep::multi_channel(
		hep::make_multi_channel_integrand<T>(
			std::ref(observables),
			generator.dimensions(),
			std::ref(generator),
			generator.map_dimensions(),
			generator.channels(),
			std::vector<hep::distribution_parameters<T>>()
		),
		std::vector<std::size_t>(1, 10)
	);

	CHECK( results.front().value() == T() );
	CHECK( results.front().error() == T() );
}

TEST_CASE("test make_mc_observables_real", "[mc_observables]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};
	std::size_t count = 2;
	T const min_energy = T(10.0);
	T const cmf_energy = T(1000.0);

	auto generator = hep::mc_phase_space_adapter<T>(
		hep::make_lusifer_phase_space_generator<T>(
			min_energy,
			cmf_energy,
			"el~el mu mu~ga ",
			constants
		)
	);

	hep::mc_observables<T> observables{hep::make_observables_real(
		test_matrix_elements<T>(set, count, 2, 4, 3),
		test_subtraction<T>(2, 4, 3),
		test_cuts<T>(count, false),
		test_recombiner<T>(count),
		test_pdf<T>(),
		test_scale_setter<T>(),
		hep::mc_trivial_distributions<T>(),
		set,
		T(1.0),
		T(0.0)
	)};

	auto const results = hep::multi_channel(
		hep::make_multi_channel_integrand<T>(
			std::ref(observables),
			generator.dimensions(),
			std::ref(generator),
			generator.map_dimensions(),
			generator.channels(),
			std::vector<hep::distribution_parameters<T>>()
		),
		std::vector<std::size_t>(1, 10)
	);

	CHECK( results.front().value() != T() );
	CHECK( results.front().error() != T() );
}
