#include "hep/ps/generate_dipole.hpp"

#include "catch.hpp"

TEST_CASE("test electroweak two jet production")
{
	std::vector<std::array<std::size_t, 3>> dipoles;

	// d d -> d d a
	std::vector<int> process = { 1, 1, 1, 1, 22 };

	// EW production at O(a^3)
	auto order = hep::coupling_order(3, 0);

	// EW dipoles
	auto type = hep::correction_type::ew;

	// all twelve dipoles where the photon is unresolved
	dipoles = {
		{ 0, 4, 1 }, { 0, 4, 2 }, { 0, 4, 3 },
		{ 1, 4, 0 }, { 1, 4, 2 }, { 1, 4, 3 },
		{ 2, 4, 0 }, { 2, 4, 1 }, { 2, 4, 3 },
		{ 3, 4, 0 }, { 3, 4, 1 }, { 3, 4, 2 }
	};

	REQUIRE( dipoles.size() == 12 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, type,
			dipole[0], dipole[1], dipole[2]);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{1, 1, 1, 1}) );
	}

	dipoles = { { 0, 2, 1 }, { 0, 2, 3 }, { 0, 2, 4 } };

	REQUIRE( dipoles.size() == 3 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, type,
			dipole[0], dipole[1], dipole[2]);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{22, 1, 1, 22}) );
	}

	dipoles = { { 1, 2, 0 }, { 1, 2, 3 }, { 1, 2, 4 } };

	REQUIRE( dipoles.size() == 3 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, type,
			dipole[0], dipole[1], dipole[2]);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{1, 22, 1, 22}) );
	}

	// non-existent dipoles
	dipoles = {
		{ 3, 2, 0 }, { 3, 2, 1 }, { 3, 2, 4 },
		{ 4, 2, 0 }, { 4, 2, 1 }, { 4, 2, 3 }
	};

	REQUIRE( dipoles.size() == 6 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, type,
			dipole[0], dipole[1], dipole[2]);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{}) );
	}

	dipoles = { { 0, 3, 1 }, { 0, 3, 2 }, { 0, 3, 4 } };

	REQUIRE( dipoles.size() == 3 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, type,
			dipole[0], dipole[1], dipole[2]);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{22, 1, 1, 22}) );
	}

	dipoles = { { 1, 3, 0 }, { 1, 3, 2 }, { 1, 3, 4 } };

	REQUIRE( dipoles.size() == 3 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, type,
			dipole[0], dipole[1], dipole[2]);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{1, 22, 1, 22}) );
	}

	// non-existent dipoles
	dipoles = {
		{ 2, 3, 0 }, { 2, 3, 1 }, { 2, 3, 4 },
		{ 4, 3, 0 }, { 4, 3, 1 }, { 4, 3, 2 }
	};

	REQUIRE( dipoles.size() == 6 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, type,
			dipole[0], dipole[1], dipole[2]);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{}) );
	}

	// there are no QCD dipoles at this order
	dipoles = {
		{ 0, 2, 1 }, { 0, 2, 3 }, { 0, 2, 4 },
		{ 1, 2, 0 }, { 1, 2, 3 }, { 1, 2, 4 },
		{ 3, 2, 0 }, { 3, 2, 1 }, { 3, 2, 4 },
		{ 4, 2, 0 }, { 4, 2, 1 }, { 4, 2, 3 },
		{ 0, 3, 1 }, { 0, 3, 2 }, { 0, 3, 4 },
		{ 1, 3, 0 }, { 1, 3, 2 }, { 1, 3, 4 },
		{ 2, 3, 0 }, { 2, 3, 1 }, { 2, 3, 4 },
		{ 4, 3, 0 }, { 4, 3, 1 }, { 4, 3, 2 },
		{ 0, 4, 1 }, { 0, 4, 2 }, { 0, 4, 3 },
		{ 1, 4, 0 }, { 1, 4, 2 }, { 1, 4, 3 },
		{ 2, 4, 0 }, { 2, 4, 1 }, { 2, 4, 3 },
		{ 3, 4, 0 }, { 3, 4, 1 }, { 3, 4, 2 }
	};

	REQUIRE( dipoles.size() == 36 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order,
			hep::correction_type::qcd, dipole[0], dipole[1], dipole[2]);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{}) );
	}
}
