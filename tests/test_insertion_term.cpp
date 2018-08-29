#include "hep/ps/correction_type.hpp"
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/particle_type.hpp"

#include <catch.hpp>

TEST_CASE("test insertion_term members", "[insertion_term]")
{
    constexpr auto fermion = hep::particle_type::fermion;

    hep::insertion_term term0{0, hep::correction_type::qcd};

    CHECK( term0.type() == hep::insertion_term_type::born );
    CHECK( term0.initial_particle() == 0 );
    CHECK( term0.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term1{0, fermion, 1, hep::correction_type::qcd};

    CHECK( term1.emitter() == 0 );
    CHECK( term1.spectator() == 1 );
    CHECK( term1.emitter_type() == fermion );
    CHECK( term1.type() == hep::insertion_term_type::initial_initial );
    CHECK( term1.initial_particle() == 0 );
    CHECK( term1.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term2{0, fermion, 2, hep::correction_type::qcd};

    CHECK( term2.emitter() == 0 );
    CHECK( term2.spectator() == 2 );
    CHECK( term2.emitter_type() == fermion );
    CHECK( term2.type() == hep::insertion_term_type::initial_final );
    CHECK( term2.initial_particle() == 0 );
    CHECK( term2.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term3{0, fermion, 3, hep::correction_type::qcd};

    CHECK( term3.emitter() == 0 );
    CHECK( term3.spectator() == 3 );
    CHECK( term3.emitter_type() == fermion );
    CHECK( term3.type() == hep::insertion_term_type::initial_final );
    CHECK( term3.initial_particle() == 0 );
    CHECK( term3.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term4{1, fermion, 0, hep::correction_type::qcd};

    CHECK( term4.emitter() == 1 );
    CHECK( term4.spectator() == 0 );
    CHECK( term4.emitter_type() == fermion );
    CHECK( term4.type() == hep::insertion_term_type::initial_initial );
    CHECK( term4.initial_particle() == 1 );
    CHECK( term4.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term5{1, fermion, 2, hep::correction_type::qcd};

    CHECK( term5.emitter() == 1 );
    CHECK( term5.spectator() == 2 );
    CHECK( term5.emitter_type() == fermion );
    CHECK( term5.type() == hep::insertion_term_type::initial_final );
    CHECK( term5.initial_particle() == 1 );
    CHECK( term5.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term6{1, fermion, 3, hep::correction_type::qcd};

    CHECK( term6.emitter() == 1 );
    CHECK( term6.spectator() == 3 );
    CHECK( term6.emitter_type() == fermion );
    CHECK( term6.type() == hep::insertion_term_type::initial_final );
    CHECK( term6.initial_particle() == 1 );
    CHECK( term6.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term7{2, fermion, 0, hep::correction_type::qcd};

    CHECK( term7.emitter() == 2 );
    CHECK( term7.spectator() == 0 );
    CHECK( term7.emitter_type() == fermion );
    CHECK( term7.type() == hep::insertion_term_type::final_initial );
    CHECK( term7.initial_particle() == 0 );
    CHECK( term7.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term8{2, fermion, 1, hep::correction_type::qcd};

    CHECK( term8.emitter() == 2 );
    CHECK( term8.spectator() == 1 );
    CHECK( term8.emitter_type() == fermion );
    CHECK( term8.type() == hep::insertion_term_type::final_initial );
    CHECK( term8.initial_particle() == 1 );
    CHECK( term8.corr_type() == hep::correction_type::qcd );

    hep::insertion_term term9{2, fermion, 3, hep::correction_type::qcd};

    CHECK( term9.emitter() == 2 );
    CHECK( term9.spectator() == 3 );
    CHECK( term9.emitter_type() == fermion );
    CHECK( term9.type() == hep::insertion_term_type::final_final );
    CHECK( term9.corr_type() == hep::correction_type::qcd );
}
