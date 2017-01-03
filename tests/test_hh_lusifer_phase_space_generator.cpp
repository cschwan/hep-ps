#include "hep/ps/hh_phase_space_generator.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"

#include "catch.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

using T = double;

hep::lusifer_constants<T> constants(
	T(125.09), T(4.0e-3),
	T(174.2), T(1.41),
	T(80.385), T(2.085),
	T(91.1876), T(2.4952)
);

TEST_CASE("constructors", "[hh_lusifer_phase_space_generator]")
{
	// e+e- -> muon pair
	hep::hh_lusifer_phase_space_generator<T> psg1(
		T(10.0),
		"el~el mu mu~",
		constants
	);

	// two s-channels with a photon or Z boson
	REQUIRE( psg1.channels()       ==  2 );
	REQUIRE( psg1.dimensions()     ==  4 );
	REQUIRE( psg1.map_dimensions() == 16 );

	hep::hh_lusifer_phase_space_generator<T> psg2(
		T(10.0),
		"sq~uq W+ W+ dq cq~",
		constants
	);

	REQUIRE( psg2.channels()       == 43 );
	REQUIRE( psg2.dimensions()     == 10 );
	REQUIRE( psg2.map_dimensions() == 24 );

	// pp -> 2 jets + two pairs of leptons
	hep::hh_lusifer_phase_space_generator<T> psg3(
		T(10.0),
		"sq~uq ne el~nm mu~dq cq~",
		constants
	);

	REQUIRE( psg3.channels()       == 93 );
	REQUIRE( psg3.dimensions()     == 16 );
	REQUIRE( psg3.map_dimensions() == 32 );
}

TEST_CASE("phase space generation", "[hh_lusifer_phase_space_generator]")
{
	hep::hh_lusifer_phase_space_generator<T> psg(
		T(10.0),
		"sq~uq ne el~nm mu~dq cq~",
		constants
	);

	std::mt19937 rng;
	std::vector<T> random_numbers(psg.dimensions());

	// energy should not be much smaller than the masses we specified
	T const energy = T(1000.0);

	std::vector<T> p(psg.map_dimensions());
	std::vector<T> densities(psg.channels());

	std::generate(random_numbers.begin(), random_numbers.end(), [&](){
		return std::generate_canonical<T, std::numeric_limits<T>::digits>(rng);
	});

	// random channel
	std::size_t const channel = 53;

	psg.generate(random_numbers, p, energy, channel);

	// check if the momentum fractions are between zero and one
	CHECK( psg.info().x1() >= T() );
	CHECK( psg.info().x1() < T(1.0) );
	CHECK( psg.info().x2() >= T() );
	CHECK( psg.info().x2() < T(1.0) );

	T const shat = (p.at(0) + p.at(0)) * (p.at(0) + p.at(0));
	T const s = shat / (psg.info().x1() * psg.info().x2());

	// check if the computed energy and momentum fractions are consistent
	CHECK( s == Approx(energy * energy) );
}
