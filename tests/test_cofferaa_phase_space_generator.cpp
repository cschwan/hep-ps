#include "hep/ps/cofferaa_phase_space_generator.hpp"

#include "catch.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

using T = double;

// TODO: reference implementation does not pass the test when massive particles
// have zero width

hep::lusifer_constants<T> constants(
	T(125.09), T(4.0e-3),
	T(174.2), T(1.41),
	T(80.385), T(2.085),
	T(91.1876), T(2.4952)
);

TEST_CASE("constructors", "[cofferaa_phase_space_generator]")
{
	T const min_energy = T(10.0);
	T const cmf_energy = T(1000.0);

	// e+e- -> muon pair
	auto psg1 = hep::make_cofferaa_phase_space_generator(
		min_energy,
		cmf_energy,
		std::vector<int>{11, -11, 13, -13},
		constants
	);

	// two s-channels with a photon or Z boson
	CHECK( psg1->channels()       ==  2 );
	CHECK( psg1->dimensions()     ==  4 );
	CHECK( psg1->map_dimensions() == 16 );

	// e+e- -> up-quark pair
	auto psg2 = hep::make_cofferaa_phase_space_generator(
		min_energy,
		cmf_energy,
		std::vector<int>{11, -11, 2, -2},
		constants
	);

	// two s-channels with a photon or Z boson
	CHECK( psg2->channels()       ==  2 );
	CHECK( psg2->dimensions()     ==  4 );
	CHECK( psg2->map_dimensions() == 16 );

	// pp -> 2 jets + two pairs of leptons
	auto psg3 = hep::make_cofferaa_phase_space_generator(
		min_energy,
		cmf_energy,
		std::vector<int>{-3, 2, 12, -11, 14, -13, 1, -4},
		constants
	);

	// TODO: why are there one channel less than what LUSIFER's PSG returns?
	CHECK( psg3->channels()       == 93-1 );
	CHECK( psg3->dimensions()     == 16 );
	CHECK( psg3->map_dimensions() == 32 );

	std::vector<std::tuple<int, int, int>> dipoles = {
		std::make_tuple(1, 9, 8),
		std::make_tuple(2, 9, 7),
		std::make_tuple(7, 9, 2),
		std::make_tuple(8, 9, 1)
	};

	// pp -> 2 jets + two pairs of leptons + gluon
	auto psg4 = hep::make_cofferaa_phase_space_generator(
		min_energy,
		cmf_energy,
		std::vector<int>{-3, 2, 12, -11, 14, -13, 1, -4, 26},
		constants,
		dipoles
	);

	CHECK( psg4->channels()       == (456+4*92) );
	CHECK( psg4->dimensions()     == 19 );
	CHECK( psg4->map_dimensions() == 36 );

	// pp -> 2 jets + two pairs of leptons
	auto psg5 = hep::make_minimal_cofferaa_phase_space_generator(
		min_energy,
		cmf_energy,
		std::vector<int>{-3, 2, 12, -11, 14, -13, 1, -4},
		constants
	);

	// TODO: why are there 16 and not 8 channels?
	CHECK( psg5->channels()       == 16 );
	CHECK( psg5->dimensions()     == 16 );
	CHECK( psg5->map_dimensions() == 32 );

	// pp -> 2 jets + two pairs of leptons + gluon
	auto psg6 = hep::make_minimal_cofferaa_phase_space_generator(
		min_energy,
		cmf_energy,
		std::vector<int>{-3, 2, 12, -11, 14, -13, 1, -4, 26},
		constants,
		dipoles
	);

	// TODO: why are there 96 and not 56/2*56 channels?
	CHECK( psg6->channels()       == (96+4*16) );
	CHECK( psg6->dimensions()     == 19 );
	CHECK( psg6->map_dimensions() == 36 );
}
