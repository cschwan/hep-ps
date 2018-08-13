#include "hep/ps/list_phase_space_generator.hpp"
#include "hep/ps/luminosity_info.hpp"

#include <catch.hpp>

#include <vector>

using T = HEP_TYPE_T;

TEST_CASE("test list_phase_space_generator", "[list_phase_space_generator]")
{
	std::vector<hep::list_phase_space_point<T>> list = {
		hep::list_phase_space_point<T>{
			hep::luminosity_info<T>{T(0.1), T(0.2), T(100.0), T(1.0)},
			std::vector<T>{ 1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0},
		},
		hep::list_phase_space_point<T>{
			hep::luminosity_info<T>{T(0.2), T(0.3), T(200.0), T(2.0)},
			std::vector<T>{ 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0},
		},
		hep::list_phase_space_point<T>{
			hep::luminosity_info<T>{T(0.3), T(0.4), T(300.0), T(3.0)},
			std::vector<T>{17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0},
		},
		hep::list_phase_space_point<T>{
			hep::luminosity_info<T>{T(0.4), T(0.5), T(400.0), T(4.0)},
			std::vector<T>{25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0}
		}
	};

	auto psg = hep::make_list_phase_space_generator(list);

	REQUIRE( psg->channels() == 1 );
	REQUIRE( psg->dimensions() == 0 );
	REQUIRE( psg->map_dimensions() == 8 );

	std::vector<T> ps;

	psg->generate(std::vector<T>{}, ps, 0);

	REQUIRE( ps.size() == 8 );

	CHECK( ps.at(0) == T( 1.0) );
	CHECK( ps.at(1) == T( 2.0) );
	CHECK( ps.at(2) == T( 3.0) );
	CHECK( ps.at(3) == T( 4.0) );
	CHECK( ps.at(4) == T( 5.0) );
	CHECK( ps.at(5) == T( 6.0) );
	CHECK( ps.at(6) == T( 7.0) );
	CHECK( ps.at(7) == T( 8.0) );

	CHECK( psg->info().x1() == T(0.1) );
	CHECK( psg->info().x2() == T(0.2) );
	CHECK( psg->info().energy_squared() == T(100.0) );
	CHECK( psg->info().rapidity_shift() == T(1.0) );

	psg->generate(std::vector<T>{}, ps, 0);

	REQUIRE( ps.size() == 8 );

	CHECK( ps.at(0) == T( 9.0) );
	CHECK( ps.at(1) == T(10.0) );
	CHECK( ps.at(2) == T(11.0) );
	CHECK( ps.at(3) == T(12.0) );
	CHECK( ps.at(4) == T(13.0) );
	CHECK( ps.at(5) == T(14.0) );
	CHECK( ps.at(6) == T(15.0) );
	CHECK( ps.at(7) == T(16.0) );

	CHECK( psg->info().x1() == T(0.2) );
	CHECK( psg->info().x2() == T(0.3) );
	CHECK( psg->info().energy_squared() == T(200.0) );
	CHECK( psg->info().rapidity_shift() == T(2.0) );

	psg->generate(std::vector<T>{}, ps, 0);

	REQUIRE( ps.size() == 8 );

	CHECK( ps.at(0) == T(17.0) );
	CHECK( ps.at(1) == T(18.0) );
	CHECK( ps.at(2) == T(19.0) );
	CHECK( ps.at(3) == T(20.0) );
	CHECK( ps.at(4) == T(21.0) );
	CHECK( ps.at(5) == T(22.0) );
	CHECK( ps.at(6) == T(23.0) );
	CHECK( ps.at(7) == T(24.0) );

	CHECK( psg->info().x1() == T(0.3) );
	CHECK( psg->info().x2() == T(0.4) );
	CHECK( psg->info().energy_squared() == T(300.0) );
	CHECK( psg->info().rapidity_shift() == T(3.0) );

	psg->generate(std::vector<T>{}, ps, 0);

	REQUIRE( ps.size() == 8 );

	CHECK( ps.at(0) == T(25.0) );
	CHECK( ps.at(1) == T(26.0) );
	CHECK( ps.at(2) == T(27.0) );
	CHECK( ps.at(3) == T(28.0) );
	CHECK( ps.at(4) == T(29.0) );
	CHECK( ps.at(5) == T(30.0) );
	CHECK( ps.at(6) == T(31.0) );

	CHECK( psg->info().x1() == T(0.4) );
	CHECK( psg->info().x2() == T(0.5) );
	CHECK( psg->info().energy_squared() == T(400.0) );
	CHECK( psg->info().rapidity_shift() == T(4.0) );

	psg->generate(std::vector<T>{}, ps, 0);

	REQUIRE( ps.size() == 8 );

	CHECK( ps.at(0) == T( 1.0) );
	CHECK( ps.at(1) == T( 2.0) );
	CHECK( ps.at(2) == T( 3.0) );
	CHECK( ps.at(3) == T( 4.0) );
	CHECK( ps.at(4) == T( 5.0) );
	CHECK( ps.at(5) == T( 6.0) );
	CHECK( ps.at(6) == T( 7.0) );
	CHECK( ps.at(7) == T( 8.0) );

	CHECK( psg->info().x1() == T(0.1) );
	CHECK( psg->info().x2() == T(0.2) );
	CHECK( psg->info().energy_squared() == T(100.0) );
	CHECK( psg->info().rapidity_shift() == T(1.0) );
}
