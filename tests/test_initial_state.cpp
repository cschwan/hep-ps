#include "hep/ps/initial_state.hpp"

#include "catch.hpp"

TEST_CASE("test partons_in_initial_state_set", "[initial_state]")
{
	hep::initial_state_set seti1{hep::initial_state::q43_uu};
	hep::parton_set setp1{hep::parton::up};

	CHECK( hep::partons_in_initial_state_set(seti1) == setp1 );

	hep::initial_state_set seti2{hep::initial_state::q43_cu};
	hep::parton_set setp2{hep::parton::up, hep::parton::charm};

	CHECK( hep::partons_in_initial_state_set(seti2) == setp2 );

	hep::initial_state_set seti3{
		hep::initial_state::q43_cu,
		hep::initial_state::q23_sd,
		hep::initial_state::q23_ug
	};
	hep::parton_set setp3{
		hep::parton::up,
		hep::parton::charm,
		hep::parton::anti_strange,
		hep::parton::anti_down,
		hep::parton::gluon
	};

	CHECK( hep::partons_in_initial_state_set(seti3) == setp3 );
}
