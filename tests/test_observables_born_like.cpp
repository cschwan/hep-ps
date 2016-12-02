#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/observables_born_like.hpp"

#include "test_structures.hpp"

#include "catch.hpp"

#include <cstddef>

using T = HEP_TYPE_T;

void test_observables_born_like(hep::initial_state_set set, std::size_t count)
{
	auto observables_born_like = hep::make_observables_born_like<T>(
		test_matrix_elements<T>(set, count, 0, 0, 0),
		test_cuts<T>(count),
		test_recombiner<T>(count),
		test_luminosities<T>(set),
		test_scale_setter<T>(),
		T(1.0)
	);

	std::vector<T> phase_space_point(4 * (count + 2));
	hep::luminosity_info<T> info{T(0.5), T(0.5), T(1024.0), T()};

	T const test_result = T(0.5) / T(1024.0);
	T const result = observables_born_like(phase_space_point, info, set);

	REQUIRE( result == test_result );
}

TEST_CASE("observables_born_like", "[observables_born_like]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};

	test_observables_born_like(set, 4);
	test_observables_born_like(set, 4);
}

