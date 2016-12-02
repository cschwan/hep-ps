#include "hep/ps/initial_state.hpp"

#include "catch.hpp"

TEST_CASE("pos or neg shift", "[initial_state]")
{
	REQUIRE( hep::state_has_pos_shift(hep::initial_state::q43_cu) == true );
	REQUIRE( hep::state_has_pos_shift(hep::initial_state::q33_du) == true );
	REQUIRE( hep::state_has_pos_shift(hep::initial_state::q33_su) == true );
	REQUIRE( hep::state_has_pos_shift(hep::initial_state::q33_dc) == true );
	REQUIRE( hep::state_has_pos_shift(hep::initial_state::q33_sc) == true );
	REQUIRE( hep::state_has_pos_shift(hep::initial_state::q23_sd) == true );

	REQUIRE( hep::state_has_neg_shift(hep::initial_state::q43_uc) == true );
	REQUIRE( hep::state_has_neg_shift(hep::initial_state::q33_ud) == true );
	REQUIRE( hep::state_has_neg_shift(hep::initial_state::q33_us) == true );
	REQUIRE( hep::state_has_neg_shift(hep::initial_state::q33_cd) == true );
	REQUIRE( hep::state_has_neg_shift(hep::initial_state::q33_cs) == true );
	REQUIRE( hep::state_has_neg_shift(hep::initial_state::q23_ds) == true );
}
