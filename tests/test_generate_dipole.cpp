#include "hep/ps/correction_type.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/generate_dipole.hpp"
#include "hep/ps/particle_type.hpp"

#include "catch.hpp"

TEST_CASE("test electroweak two jet production")
{
	std::vector<hep::dipole> dipoles;

	// d d -> d d a
	std::vector<int> process = { 1, 1, 1, 1, 22 };

	// EW production at O(a^3)
	auto order = hep::coupling_order(3, 0);

	constexpr auto f = hep::particle_type::fermion;
	constexpr auto b = hep::particle_type::boson;
	constexpr auto ew = hep::correction_type::ew;
	constexpr auto qcd = hep::correction_type::qcd;

	// all twelve dipoles where the photon is unresolved
	dipoles = {
		hep::dipole{0, 4, 1, f, b, f, ew}, hep::dipole{0, 4, 2, f, b, f, ew},
		hep::dipole{0, 4, 3, f, b, f, ew}, hep::dipole{1, 4, 0, f, b, f, ew},
		hep::dipole{1, 4, 2, f, b, f, ew}, hep::dipole{1, 4, 3, f, b, f, ew},
		hep::dipole{2, 4, 0, f, b, f, ew}, hep::dipole{2, 4, 1, f, b, f, ew},
		hep::dipole{2, 4, 3, f, b, f, ew}, hep::dipole{3, 4, 0, f, b, f, ew},
		hep::dipole{3, 4, 1, f, b, f, ew}, hep::dipole{3, 4, 2, f, b, f, ew}
	};

	REQUIRE( dipoles.size() == 12 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, dipole);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{1, 1, 1, 1}) );
	}

	dipoles = {
		hep::dipole{0, 2, 1, f, f, f, ew},
		hep::dipole{0, 2, 3, f, f, f, ew},
		hep::dipole{0, 2, 4, f, f, b, ew}
	};

	REQUIRE( dipoles.size() == 3 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, dipole);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{22, 1, 1, 22}) );
	}

	dipoles = {
		hep::dipole{1, 2, 0, f, f, f, ew},
		hep::dipole{1, 2, 3, f, f, f, ew},
		hep::dipole{1, 2, 4, f, f, b, ew}
	};

	REQUIRE( dipoles.size() == 3 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, dipole);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{1, 22, 1, 22}) );
	}

	// non-existent dipoles
	dipoles = {
		hep::dipole{3, 2, 0, f, f, f, ew}, hep::dipole{3, 2, 1, f, f, f, ew},
		hep::dipole{3, 2, 4, f, f, b, ew}, hep::dipole{4, 2, 0, b, f, f, ew},
		hep::dipole{4, 2, 1, b, f, f, ew}, hep::dipole{4, 2, 3, b, f, f, ew}
	};

	REQUIRE( dipoles.size() == 6 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, dipole);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{}) );
	}

	dipoles = {
		hep::dipole{0, 3, 1, f, f, f, ew},
		hep::dipole{0, 3, 2, f, f, f, ew},
		hep::dipole{0, 3, 4, f, f, b, ew}
	};

	REQUIRE( dipoles.size() == 3 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, dipole);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{22, 1, 1, 22}) );
	}

	dipoles = {
		hep::dipole{1, 3, 0, f, f, f, ew},
		hep::dipole{1, 3, 2, f, f, f, ew},
		hep::dipole{1, 3, 4, f, f, b, ew}
	};

	REQUIRE( dipoles.size() == 3 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, dipole);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{1, 22, 1, 22}) );
	}

	// non-existent dipoles
	dipoles = {
		hep::dipole{2, 3, 0, f, f, f, ew}, hep::dipole{2, 3, 1, f, f, f, ew},
		hep::dipole{2, 3, 4, f, f, b, ew}, hep::dipole{4, 3, 0, b, f, f, ew},
		hep::dipole{4, 3, 1, b, f, f, ew}, hep::dipole{4, 3, 2, b, f, f, ew}
	};

	REQUIRE( dipoles.size() == 6 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, dipole);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{}) );
	}

	// there are no QCD dipoles at this order
	dipoles = {
		hep::dipole{0, 2, 1, f, f, f, qcd}, {0, 2, 3, f, f, f, qcd},
		hep::dipole{0, 2, 4, f, f, b, qcd}, {1, 2, 0, f, f, f, qcd},
		hep::dipole{1, 2, 3, f, f, f, qcd}, {1, 2, 4, f, f, b, qcd},
		hep::dipole{3, 2, 0, f, f, f, qcd}, {3, 2, 1, f, f, f, qcd},
		hep::dipole{3, 2, 4, f, f, b, qcd}, {4, 2, 0, b, f, f, qcd},
		hep::dipole{4, 2, 1, b, f, f, qcd}, {4, 2, 3, b, f, f, qcd},
		hep::dipole{0, 3, 1, f, f, f, qcd}, {0, 3, 2, f, f, f, qcd},
		hep::dipole{0, 3, 4, f, f, b, qcd}, {1, 3, 0, f, f, f, qcd},
		hep::dipole{1, 3, 2, f, f, f, qcd}, {1, 3, 4, f, f, b, qcd},
		hep::dipole{2, 3, 0, f, f, f, qcd}, {2, 3, 1, f, f, f, qcd},
		hep::dipole{2, 3, 4, f, f, b, qcd}, {4, 3, 0, b, f, f, qcd},
		hep::dipole{4, 3, 1, b, f, f, qcd}, {4, 3, 2, b, f, f, qcd},
		hep::dipole{0, 4, 1, f, b, f, qcd}, {0, 4, 2, f, b, f, qcd},
		hep::dipole{0, 4, 3, f, b, f, qcd}, {1, 4, 0, f, b, f, qcd},
		hep::dipole{1, 4, 2, f, f, f, qcd}, {1, 4, 3, f, b, f, qcd},
		hep::dipole{2, 4, 0, f, b, f, qcd}, {2, 4, 1, f, b, f, qcd},
		hep::dipole{2, 4, 3, f, b, f, qcd}, {3, 4, 0, f, b, f, qcd},
		hep::dipole{3, 4, 1, f, b, f, qcd}, {3, 4, 2, f, b, f, qcd}
	};

	REQUIRE( dipoles.size() == 36 );

	for (auto const& dipole : dipoles)
	{
		auto const& result = hep::generate_dipole(process, order, dipole);

		CHECK_THAT( result , Catch::Equals(std::vector<int>{}) );
	}
}
