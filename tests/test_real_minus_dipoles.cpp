#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/real_minus_dipoles.hpp"

#include "test_structures.hpp"

#include "catch.hpp"

#include <cstddef>

using T = HEP_TYPE_T;

void test_real_minus_dipoles(
	hep::initial_state_set set,
	std::size_t count,
	std::size_t emitter,
	std::size_t unresolved,
	std::size_t spectator,
	bool inclusive
) {
	auto real_minus_dipoles = hep::make_real_minus_dipoles<T>(
		test_matrix_elements<T>(set, count, emitter, unresolved, spectator),
		test_subtraction<T>(emitter, unresolved, spectator),
		test_cuts<T>(count),
		test_recombiner<T>(count),
		inclusive
	);

	std::vector<T> phase_space(4 * (count + 2));

	T const test_result = T(inclusive ? 0.5 : -0.5);
	auto result = real_minus_dipoles(phase_space, T(), set);

	for (auto const process : hep::initial_state_list())
	{
		if (set.includes(process))
		{
			REQUIRE( result.get(process) == test_result );
		}
		else
		{
			REQUIRE( result.get(process) == T() );
		}
	}
}

TEST_CASE("real_minus_dipoles", "[real_minus_dipoles]")
{
	hep::initial_state_set set;
	set.add(hep::initial_state::q43_cu);

	test_real_minus_dipoles(set, 4, 2, 3, 0, true);
	test_real_minus_dipoles(set, 4, 2, 3, 0, false);
}
