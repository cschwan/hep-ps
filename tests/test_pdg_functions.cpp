#include "hep/ps/pdg_functions.hpp"

#include "catch.hpp"

TEST_CASE("test parton_to_pdg_id()", "[pdg_functions]")
{
	// http://pdg.lbl.gov/2017/reviews/rpp2017-rev-monte-carlo-numbering.pdf
	CHECK( hep::parton_to_pdg_id(hep::parton::anti_charm)   == -4 );
	CHECK( hep::parton_to_pdg_id(hep::parton::anti_strange) == -3 );
	CHECK( hep::parton_to_pdg_id(hep::parton::anti_up)      == -2 );
	CHECK( hep::parton_to_pdg_id(hep::parton::anti_down)    == -1 );
	CHECK( hep::parton_to_pdg_id(hep::parton::down)         ==  1 );
	CHECK( hep::parton_to_pdg_id(hep::parton::up)           ==  2 );
	CHECK( hep::parton_to_pdg_id(hep::parton::strange)      ==  3 );
	CHECK( hep::parton_to_pdg_id(hep::parton::charm)        ==  4 );
	CHECK( hep::parton_to_pdg_id(hep::parton::gluon)        == 21 );
	CHECK( hep::parton_to_pdg_id(hep::parton::photon)       == 22 );
}

TEST_CASE("test ", "[pdg_functions]")
{
	std::string const input = "-2 -4 -> -1 -3 11 -12 13 -14 22";
	std::string const output = hep::pdg_ids_to_ol_process_string(
		hep::ol_process_string_to_pdg_ids(input));

	CHECK( input == output );
}
