#include "hep/ps/insertion_term.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/particle_type.hpp"

#include "catch.hpp"

TEST_CASE("test insertion_term members", "[insertion_term]")
{
	constexpr auto boson = hep::particle_type::boson;
	constexpr auto fermion = hep::particle_type::fermion;

	hep::insertion_term term0;

	CHECK( term0.type() == hep::insertion_term_type::born );

	hep::insertion_term term1{0, fermion, 1};

	CHECK( term1.emitter() == 0 );
	CHECK( term1.spectator() == 1 );
	CHECK( term1.emitter_type() == fermion );
	CHECK( term1.type() == hep::insertion_term_type::initial_initial );

	hep::insertion_term term2{0, fermion, 2};

	CHECK( term2.emitter() == 0 );
	CHECK( term2.spectator() == 2 );
	CHECK( term2.emitter_type() == fermion );
	CHECK( term2.type() == hep::insertion_term_type::initial_final );

	hep::insertion_term term3{0, fermion, 3};

	CHECK( term3.emitter() == 0 );
	CHECK( term3.spectator() == 3 );
	CHECK( term3.emitter_type() == fermion );
	CHECK( term3.type() == hep::insertion_term_type::initial_final );

	hep::insertion_term term4{1, fermion, 0};

	CHECK( term4.emitter() == 1 );
	CHECK( term4.spectator() == 0 );
	CHECK( term4.emitter_type() == fermion );
	CHECK( term4.type() == hep::insertion_term_type::initial_initial );

	hep::insertion_term term5{1, fermion, 2};

	CHECK( term5.emitter() == 1 );
	CHECK( term5.spectator() == 2 );
	CHECK( term5.emitter_type() == fermion );
	CHECK( term5.type() == hep::insertion_term_type::initial_final );

	hep::insertion_term term6{1, fermion, 3};

	CHECK( term6.emitter() == 1 );
	CHECK( term6.spectator() == 3 );
	CHECK( term6.emitter_type() == fermion );
	CHECK( term6.type() == hep::insertion_term_type::initial_final );

	hep::insertion_term term7{2, fermion, 0};

	CHECK( term7.emitter() == 2 );
	CHECK( term7.spectator() == 0 );
	CHECK( term7.emitter_type() == fermion );
	CHECK( term7.type() == hep::insertion_term_type::final_initial );

	hep::insertion_term term8{2, fermion, 1};

	CHECK( term8.emitter() == 2 );
	CHECK( term8.spectator() == 1 );
	CHECK( term8.emitter_type() == fermion );
	CHECK( term8.type() == hep::insertion_term_type::final_initial );

	hep::insertion_term term9{2, fermion, 3};

	CHECK( term9.emitter() == 2 );
	CHECK( term9.spectator() == 3 );
	CHECK( term9.emitter_type() == fermion );
	CHECK( term9.type() == hep::insertion_term_type::final_final );
}
