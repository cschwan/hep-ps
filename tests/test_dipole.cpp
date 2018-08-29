#include "hep/ps/correction_type.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_type.hpp"
#include "hep/ps/particle_type.hpp"

#include <catch.hpp>

TEST_CASE("test dipole members", "[dipole]")
{
    constexpr auto boson = hep::particle_type::boson;
    constexpr auto fermion = hep::particle_type::fermion;

    hep::dipole dipole1{0, 4, 1, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole1.emitter() == 0 );
    CHECK( dipole1.unresolved() == 4 );
    CHECK( dipole1.spectator() == 1 );
    CHECK( dipole1.emitter_type() == fermion );
    CHECK( dipole1.unresolved_type() == boson );
    CHECK( dipole1.spectator_type() == fermion );
    CHECK( dipole1.type() == hep::dipole_type::initial_initial );
    CHECK( dipole1.corr_type() == hep::correction_type::ew );

    hep::dipole dipole2{0, 4, 2, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole2.emitter() == 0 );
    CHECK( dipole2.unresolved() == 4 );
    CHECK( dipole2.spectator() == 2 );
    CHECK( dipole2.emitter_type() == fermion );
    CHECK( dipole2.unresolved_type() == boson );
    CHECK( dipole2.spectator_type() == fermion );
    CHECK( dipole2.type() == hep::dipole_type::initial_final );
    CHECK( dipole2.corr_type() == hep::correction_type::ew );

    hep::dipole dipole3{0, 4, 3, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole3.emitter() == 0 );
    CHECK( dipole3.unresolved() == 4 );
    CHECK( dipole3.spectator() == 3 );
    CHECK( dipole3.emitter_type() == fermion );
    CHECK( dipole3.unresolved_type() == boson );
    CHECK( dipole3.spectator_type() == fermion );
    CHECK( dipole3.type() == hep::dipole_type::initial_final );
    CHECK( dipole3.corr_type() == hep::correction_type::ew );

    hep::dipole dipole4{1, 4, 0, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole4.emitter() == 1 );
    CHECK( dipole4.unresolved() == 4 );
    CHECK( dipole4.spectator() == 0 );
    CHECK( dipole4.emitter_type() == fermion );
    CHECK( dipole4.unresolved_type() == boson );
    CHECK( dipole4.spectator_type() == fermion );
    CHECK( dipole4.type() == hep::dipole_type::initial_initial );
    CHECK( dipole4.corr_type() == hep::correction_type::ew );

    hep::dipole dipole5{1, 4, 2, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole5.emitter() == 1 );
    CHECK( dipole5.unresolved() == 4 );
    CHECK( dipole5.spectator() == 2 );
    CHECK( dipole5.emitter_type() == fermion );
    CHECK( dipole5.unresolved_type() == boson );
    CHECK( dipole5.spectator_type() == fermion );
    CHECK( dipole5.type() == hep::dipole_type::initial_final );
    CHECK( dipole5.corr_type() == hep::correction_type::ew );

    hep::dipole dipole6{1, 4, 3, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole6.emitter() == 1 );
    CHECK( dipole6.unresolved() == 4 );
    CHECK( dipole6.spectator() == 3 );
    CHECK( dipole6.emitter_type() == fermion );
    CHECK( dipole6.unresolved_type() == boson );
    CHECK( dipole6.spectator_type() == fermion );
    CHECK( dipole6.type() == hep::dipole_type::initial_final );
    CHECK( dipole6.corr_type() == hep::correction_type::ew );

    hep::dipole dipole7{2, 4, 0, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole7.emitter() == 2 );
    CHECK( dipole7.unresolved() == 4 );
    CHECK( dipole7.spectator() == 0 );
    CHECK( dipole7.emitter_type() == fermion );
    CHECK( dipole7.unresolved_type() == boson );
    CHECK( dipole7.spectator_type() == fermion );
    CHECK( dipole7.type() == hep::dipole_type::final_initial );
    CHECK( dipole7.corr_type() == hep::correction_type::ew );

    hep::dipole dipole8{2, 4, 1, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole8.emitter() == 2 );
    CHECK( dipole8.unresolved() == 4 );
    CHECK( dipole8.spectator() == 1 );
    CHECK( dipole8.emitter_type() == fermion );
    CHECK( dipole8.unresolved_type() == boson );
    CHECK( dipole8.spectator_type() == fermion );
    CHECK( dipole8.type() == hep::dipole_type::final_initial );
    CHECK( dipole8.corr_type() == hep::correction_type::ew );

    hep::dipole dipole9{2, 4, 3, fermion, boson, fermion,
        hep::correction_type::ew};

    CHECK( dipole9.emitter() == 2 );
    CHECK( dipole9.unresolved() == 4 );
    CHECK( dipole9.spectator() == 3 );
    CHECK( dipole9.emitter_type() == fermion );
    CHECK( dipole9.unresolved_type() == boson );
    CHECK( dipole9.spectator_type() == fermion );
    CHECK( dipole9.type() == hep::dipole_type::final_final );
    CHECK( dipole9.corr_type() == hep::correction_type::ew );
}
