#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/observables_born.hpp"
#include "hep/ps/trivial_distributions.hpp"

#include "test_structures.hpp"

#include "catch.hpp"

#include <cstddef>
#include <vector>

using T = HEP_TYPE_T;

void test_observables_born(hep::initial_state_set set, std::size_t count)
{
	auto observables = hep::make_observables_born<T>(
		test_matrix_elements<T>(set, count, 0, 0, 0),
		test_cuts<T>(count, false),
		test_recombiner<T>(count),
		test_pdf<T>(),
		test_scale_setter<T>(),
		hep::trivial_distributions<T>(),
		T(1.0)
	);

	std::vector<T> phase_space_point(4 * (count + 2));
	hep::luminosity_info<T> info{T(0.5), T(0.5), T(1024.0), T()};

	T const test_result = T(1.0) / T(1024.0);
	T const result = observables->eval(phase_space_point, info, set);

	REQUIRE( result == test_result );
}

TEST_CASE("observables_born_like", "[observables_born_like]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};

	test_observables_born(set, 4);
	test_observables_born(set, 4);
}
