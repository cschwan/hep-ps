#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/observables_born_like.hpp"
#include "hep/ps/trivial_cutter.hpp"
#include "hep/ps/trivial_recombiner.hpp"

#include "test_structures.hpp"

#include "catch.hpp"

#include <cstddef>

using T = HEP_TYPE_T;

void test_trivial_cutter_and_recombiner(
	hep::initial_state_set set,
	std::size_t count
) {
	auto observables_born_like = hep::make_observables_born_like<T>(
		test_matrix_elements<T>(set, count, 0, 0, 0),
		hep::trivial_cutter<T>(),
		hep::trivial_recombiner<T>(),
		test_pdf<T>(),
		test_scale_setter<T>(),
		T(1.0)
	);

	std::vector<T> phase_space_point(4 * (count + 2));
	hep::luminosity_info<T> info{T(0.5), T(0.5), T(1024.0), T()};

	T const test_result = T(0.5) / T(1024.0);
	T const result = observables_born_like(phase_space_point, info, set);

	REQUIRE( result == test_result );
}

TEST_CASE("check trivial cutter and recombiner",
	"[trivial_cutter],[trivial_recombiner]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};

	test_trivial_cutter_and_recombiner(set, 4);
}
