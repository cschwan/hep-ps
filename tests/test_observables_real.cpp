#include "hep/ps/initial_state.hpp"
#include "hep/ps/observables_real.hpp"
#include "hep/ps/trivial_distributions.hpp"

#include "test_structures.hpp"

#include "catch.hpp"

#include <cstddef>
#include <vector>

using T = HEP_TYPE_T;

void test_observables_real(
	hep::initial_state_set set,
	std::size_t count,
	std::size_t emitter,
	std::size_t unresolved,
	std::size_t spectator,
	bool inclusive
) {
	auto observables = hep::make_observables_real<T>(
		test_matrix_elements<T>(set, count, emitter, unresolved, spectator),
		test_subtraction<T>(emitter, unresolved, spectator),
		test_cuts<T>(count, inclusive),
		test_recombiner<T>(count),
		test_pdf<T>(),
		test_scale_setter<T>(),
		hep::trivial_distributions<T>(),
		set,
		T(1.0),
		T()
	);

	std::vector<T> real_phase_space_point(4 * (count + 3));
	hep::luminosity_info<T> info{T(0.5), T(0.5), T(1024.0), T()};

	T const test_result = T(inclusive ? 0.5 : -0.5) / T(1024.0);
	T const result = observables->eval(real_phase_space_point, info);

	REQUIRE( result == test_result );
}

TEST_CASE("observables_real", "[observables_real]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};

	test_observables_real(set, 4, 2, 3, 0, true);
	test_observables_real(set, 4, 2, 3, 0, false);
}
