#include "hep/ps/lusifer_phase_space_generator.hpp"

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

TEST_CASE("constructors", "[lusifer_phase_space_generator]")
{
	// e+e- -> muon pair
	hep::lusifer_phase_space_generator<T> psg1(
		"el~el mu mu~",
		constants
	);

	// two s-channels with a photon or Z boson
	REQUIRE( psg1.channels()       ==  2 );
	REQUIRE( psg1.dimensions()     ==  2 );
	REQUIRE( psg1.map_dimensions() == 16 );

	hep::lusifer_phase_space_generator<T> psg2(
		"sq~uq W+ W+ dq cq~",
		constants
	);

	REQUIRE( psg2.channels()       == 43 );
	REQUIRE( psg2.dimensions()     ==  8 );
	REQUIRE( psg2.map_dimensions() == 24 );

	// pp -> 2 jets + two pairs of leptons
	hep::lusifer_phase_space_generator<T> psg3(
		"sq~uq ne el~nm mu~dq cq~",
		constants
	);

	REQUIRE( psg3.channels()       == 93 );
	REQUIRE( psg3.dimensions()     == 14 );
	REQUIRE( psg3.map_dimensions() == 32 );
}

TEST_CASE("phase space generation", "[lusifer_phase_space_generator]")
{
	hep::lusifer_phase_space_generator<T> psg(
		"sq~uq ne el~nm mu~dq cq~",
		constants
	);

	std::mt19937 rng;
	std::vector<T> random_numbers(psg.dimensions());

	// energy should not be much smaller than the masses we specified
	T const energy = T(1000.0);

	std::vector<T> p(psg.map_dimensions());
	std::vector<T> densities(psg.channels());

	for (std::size_t i = 0; i != 1000; ++i)
	{
		std::generate(random_numbers.begin(), random_numbers.end(), [&](){
			return std::generate_canonical<T,
				std::numeric_limits<T>::digits>(rng);
		});

		for (std::size_t channel = 0; channel != psg.channels(); ++channel)
		{
			psg.generate(random_numbers, p, energy, channel);

			std::array<T, 4> sums = { T(), T(), T(), T() };

			for (std::size_t particle = 0; particle != 8; ++particle)
			{
				// count incoming particles as negative
				T const sign = (particle < 2) ? T(-1.0) : T(1.0);

				sums.at(0) += sign * p.at(4 * particle + 0);
				sums.at(1) += sign * p.at(4 * particle + 1);
				sums.at(2) += sign * p.at(4 * particle + 2);
				sums.at(3) += sign * p.at(4 * particle + 3);
			}

			CAPTURE( i );
			CAPTURE( channel );

			// check momentum conservation
			CHECK( sums.at(0) == Approx(T()) );
			CHECK( sums.at(1) == Approx(T()) );
			CHECK( sums.at(2) == Approx(T()) );
			CHECK( sums.at(3) == Approx(T()) );

			for (std::size_t particle = 0; particle != 8; ++particle)
			{
				T const invariant =
					p.at(4 * particle + 0) * p.at(4 * particle + 0) -
					p.at(4 * particle + 1) * p.at(4 * particle + 1) -
					p.at(4 * particle + 2) * p.at(4 * particle + 2) -
					p.at(4 * particle + 3) * p.at(4 * particle + 3);

				CAPTURE( i );
				CAPTURE( channel );
				CAPTURE( particle );

				// check on-shellness
				CHECK( invariant == Approx(T()) );
			}

			std::fill(densities.begin(), densities.end(), T());
			psg.densities(densities);

			for (std::size_t index = 0; index != densities.size(); ++index)
			{
				CAPTURE( i );
				CAPTURE( channel );
				CAPTURE( index );

				// check if the densities are non-zero
				CHECK( densities.at(index) != T() );
			}
		}
	}
}