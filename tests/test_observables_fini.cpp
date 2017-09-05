#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"

#include "hep/ps/initial_state.hpp"
#include "hep/ps/observables_fini.hpp"
#include "hep/ps/rambo_phase_space_generator.hpp"
#include "hep/ps/trivial_distributions.hpp"

#include "test_structures.hpp"

#include "catch.hpp"

#include <cstddef>
#include <functional>
#include <vector>

using T = HEP_TYPE_T;

void test_observables_fini(hep::initial_state_set set, std::size_t count)
{
	auto generator = hep::make_rambo_phase_space_generator<T>(
		T(1.0),
		T(100.0),
		count,
		1
	);

	auto integrand = hep::make_observables_fini<T>(
		test_matrix_elements<T>(set, count, 0, 0, 0),
		test_subtraction<T>(0, 0, 0),
		test_cuts<T>(count, false),
		test_recombiner<T>(count),
		test_pdf<T>(),
		test_scale_setter<T>(),
		hep::trivial_distributions<T>(),
		set,
		T(1.0)
	);

	auto const result = hep::multi_channel(
		hep::make_multi_channel_integrand<T>(
			std::ref(*integrand),
			generator->dimensions(),
			std::ref(*generator),
			generator->map_dimensions(),
			generator->channels(),
			test_distribution_parameters<T>()
		),
		std::vector<std::size_t>{1}
	);

	CHECK( result.front().value() == T() );
}

TEST_CASE("observables_fini", "[observables_fini]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};

	test_observables_fini(set, 4);
	test_observables_fini(set, 4);
}
