#include "hep/ps/initial_state.hpp"

#include <catch.hpp>

#include <iterator>

TEST_CASE("test partons_in_initial_state_set", "[initial_state]")
{
	// initial states are alphabetically sorted
	hep::initial_state_set set{
		hep::initial_state::cq_cq,
		hep::initial_state::cq_cx,
		hep::initial_state::cq_dq,
		hep::initial_state::cq_dx,
		hep::initial_state::cq_gl,
		hep::initial_state::cq_ph,
		hep::initial_state::cq_sq,
		hep::initial_state::cq_sx,
		hep::initial_state::cq_uq,
		hep::initial_state::cq_ux,
		hep::initial_state::cx_cq,
		hep::initial_state::cx_dx,
		hep::initial_state::cx_sx,
		hep::initial_state::cx_uq,
		hep::initial_state::dq_cq,
		hep::initial_state::dq_dx,
		hep::initial_state::dq_sx,
		hep::initial_state::dq_uq,
		hep::initial_state::dx_cq,
		hep::initial_state::dx_cx,
		hep::initial_state::dx_dq,
		hep::initial_state::dx_dx,
		hep::initial_state::dx_gl,
		hep::initial_state::dx_sq,
		hep::initial_state::dx_sx,
		hep::initial_state::dx_uq,
		hep::initial_state::dx_ux,
		hep::initial_state::gl_cq,
		hep::initial_state::gl_dx,
		hep::initial_state::gl_gl,
		hep::initial_state::gl_sx,
		hep::initial_state::gl_uq,
		hep::initial_state::ph_cq,
		hep::initial_state::ph_ph,
		hep::initial_state::ph_uq,
		hep::initial_state::sq_cq,
		hep::initial_state::sq_dx,
		hep::initial_state::sq_sx,
		hep::initial_state::sq_uq,
		hep::initial_state::sx_cq,
		hep::initial_state::sx_cx,
		hep::initial_state::sx_dq,
		hep::initial_state::sx_dx,
		hep::initial_state::sx_gl,
		hep::initial_state::sx_sq,
		hep::initial_state::sx_sx,
		hep::initial_state::sx_uq,
		hep::initial_state::sx_ux,
		hep::initial_state::uq_cq,
		hep::initial_state::uq_cx,
		hep::initial_state::uq_dq,
		hep::initial_state::uq_dx,
		hep::initial_state::uq_gl,
		hep::initial_state::uq_ph,
		hep::initial_state::uq_sq,
		hep::initial_state::uq_sx,
		hep::initial_state::uq_uq,
		hep::initial_state::uq_ux,
		hep::initial_state::ux_cq,
		hep::initial_state::ux_dx,
		hep::initial_state::ux_sx,
		hep::initial_state::ux_uq,
	};

	CHECK( *std::next(set.begin(),  0) == hep::initial_state::cq_cq );
	CHECK( *std::next(set.begin(),  1) == hep::initial_state::cq_cx );
	CHECK( *std::next(set.begin(),  2) == hep::initial_state::cq_dq );
	CHECK( *std::next(set.begin(),  3) == hep::initial_state::cq_dx );
	CHECK( *std::next(set.begin(),  4) == hep::initial_state::cq_gl );
	CHECK( *std::next(set.begin(),  5) == hep::initial_state::cq_ph );
	CHECK( *std::next(set.begin(),  6) == hep::initial_state::cq_sq );
	CHECK( *std::next(set.begin(),  7) == hep::initial_state::cq_sx );
	CHECK( *std::next(set.begin(),  8) == hep::initial_state::cq_uq );
	CHECK( *std::next(set.begin(),  9) == hep::initial_state::cq_ux );
	CHECK( *std::next(set.begin(), 10) == hep::initial_state::cx_cq );
	CHECK( *std::next(set.begin(), 11) == hep::initial_state::cx_dx );
	CHECK( *std::next(set.begin(), 12) == hep::initial_state::cx_sx );
	CHECK( *std::next(set.begin(), 13) == hep::initial_state::cx_uq );
	CHECK( *std::next(set.begin(), 14) == hep::initial_state::dq_cq );
	CHECK( *std::next(set.begin(), 15) == hep::initial_state::dq_dx );
	CHECK( *std::next(set.begin(), 16) == hep::initial_state::dq_sx );
	CHECK( *std::next(set.begin(), 17) == hep::initial_state::dq_uq );
	CHECK( *std::next(set.begin(), 18) == hep::initial_state::dx_cq );
	CHECK( *std::next(set.begin(), 19) == hep::initial_state::dx_cx );
	CHECK( *std::next(set.begin(), 20) == hep::initial_state::dx_dq );
	CHECK( *std::next(set.begin(), 21) == hep::initial_state::dx_dx );
	CHECK( *std::next(set.begin(), 22) == hep::initial_state::dx_gl );
	CHECK( *std::next(set.begin(), 23) == hep::initial_state::dx_sq );
	CHECK( *std::next(set.begin(), 24) == hep::initial_state::dx_sx );
	CHECK( *std::next(set.begin(), 25) == hep::initial_state::dx_uq );
	CHECK( *std::next(set.begin(), 26) == hep::initial_state::dx_ux );
	CHECK( *std::next(set.begin(), 27) == hep::initial_state::gl_cq );
	CHECK( *std::next(set.begin(), 28) == hep::initial_state::gl_dx );
	CHECK( *std::next(set.begin(), 29) == hep::initial_state::gl_gl );
	CHECK( *std::next(set.begin(), 30) == hep::initial_state::gl_sx );
	CHECK( *std::next(set.begin(), 31) == hep::initial_state::gl_uq );
	CHECK( *std::next(set.begin(), 32) == hep::initial_state::ph_cq );
	CHECK( *std::next(set.begin(), 33) == hep::initial_state::ph_ph );
	CHECK( *std::next(set.begin(), 34) == hep::initial_state::ph_uq );
	CHECK( *std::next(set.begin(), 35) == hep::initial_state::sq_cq );
	CHECK( *std::next(set.begin(), 36) == hep::initial_state::sq_dx );
	CHECK( *std::next(set.begin(), 37) == hep::initial_state::sq_sx );
	CHECK( *std::next(set.begin(), 38) == hep::initial_state::sq_uq );
	CHECK( *std::next(set.begin(), 39) == hep::initial_state::sx_cq );
	CHECK( *std::next(set.begin(), 40) == hep::initial_state::sx_cx );
	CHECK( *std::next(set.begin(), 41) == hep::initial_state::sx_dq );
	CHECK( *std::next(set.begin(), 42) == hep::initial_state::sx_dx );
	CHECK( *std::next(set.begin(), 43) == hep::initial_state::sx_gl );
	CHECK( *std::next(set.begin(), 44) == hep::initial_state::sx_sq );
	CHECK( *std::next(set.begin(), 45) == hep::initial_state::sx_sx );
	CHECK( *std::next(set.begin(), 46) == hep::initial_state::sx_uq );
	CHECK( *std::next(set.begin(), 47) == hep::initial_state::sx_ux );
	CHECK( *std::next(set.begin(), 48) == hep::initial_state::uq_cq );
	CHECK( *std::next(set.begin(), 49) == hep::initial_state::uq_cx );
	CHECK( *std::next(set.begin(), 50) == hep::initial_state::uq_dq );
	CHECK( *std::next(set.begin(), 51) == hep::initial_state::uq_dx );
	CHECK( *std::next(set.begin(), 52) == hep::initial_state::uq_gl );
	CHECK( *std::next(set.begin(), 53) == hep::initial_state::uq_ph );
	CHECK( *std::next(set.begin(), 54) == hep::initial_state::uq_sq );
	CHECK( *std::next(set.begin(), 55) == hep::initial_state::uq_sx );
	CHECK( *std::next(set.begin(), 56) == hep::initial_state::uq_uq );
	CHECK( *std::next(set.begin(), 57) == hep::initial_state::uq_ux );
	CHECK( *std::next(set.begin(), 58) == hep::initial_state::ux_cq );
	CHECK( *std::next(set.begin(), 59) == hep::initial_state::ux_dx );
	CHECK( *std::next(set.begin(), 60) == hep::initial_state::ux_sx );
	CHECK( *std::next(set.begin(), 61) == hep::initial_state::ux_uq );
	CHECK(  std::next(set.begin(), 62) == set.end() );
}

TEST_CASE("are `state_parton_{one,two}` complete?")
{
	for (auto state : hep::initial_state_list())
	{
		CAPTURE( state );

		CHECK_NOTHROW( state_parton_one(state) );
		CHECK_NOTHROW( state_parton_two(state) );
	}
}
