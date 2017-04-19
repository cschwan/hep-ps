#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/observables_fini.hpp"
#include "hep/ps/trivial_distributions.hpp"

#include "test_structures.hpp"

#include "catch.hpp"

#include <cstddef>
#include <vector>

using T = HEP_TYPE_T;

void test_observables_fint(
	hep::initial_state_set set,
	std::size_t count
) {
	auto observables = hep::make_observables_fini<T>(
		test_matrix_elements<T>(set, count, 0, 0, 0),
		test_subtraction<T>(0, 0, 0),
		test_cuts<T>(count, false),
		test_recombiner<T>(count),
		test_pdf<T>(),
		test_scale_setter<T>(),
		hep::trivial_distributions<T>(),
		T(1.0)
	);

	std::vector<T> phase_space_point(4 * (count + 2));
	hep::luminosity_info<T> info{T(0.5), T(0.5), T(1024.0), T()};

	T const test_result = T();

	T const x = T(0.5);
	std::vector<T> numbers = { x };
	hep::random_numbers<T> rans(numbers);
	T const result = observables->eval(phase_space_point, info, rans, set);

	REQUIRE( result == test_result );
}

TEST_CASE("observables_finite_insertion", "[observables_finite_insertion]")
{
	hep::initial_state_set set{hep::initial_state::q43_cu};

	test_observables_fint(set, 4);
	test_observables_fint(set, 4);
}
