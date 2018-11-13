#include "hep/ps/lusifer_phase_space_generator.hpp"

#include <catch.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

using T = double;

// TODO: reference implementation does not pass the test when massive particles
// have zero width

hep::lusifer_constants<T> constants(
    T(125.09), T(4.0e-3),
    T(174.2), T(1.41),
    T(80.385), T(2.085),
    T(91.1876), T(2.4952)
);

TEST_CASE("constructors", "[lusifer_phase_space_generator]")
{
    T const min_energy = T(10.0);
    T const cmf_energy = T(1000.0);

    // e+e- -> muon pair
    auto psg1 = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "el~el mu mu~",
        constants
    );

    // two s-channels with a photon or Z boson
    CHECK( psg1->channels()       ==  2 );
    CHECK( psg1->dimensions()     ==  4 );
    CHECK( psg1->map_dimensions() == 16 );

    for (std::size_t extra = 0; extra != 4; ++extra)
    {
        // e+e- -> muon pair + extra random numbers
        auto psg = hep::make_lusifer_phase_space_generator(
            min_energy,
            cmf_energy,
            "el~el mu mu~",
            constants,
            extra
        );

        CHECK( psg->channels()       ==  2 );
        CHECK( psg->dimensions()     ==  4 + extra );
        CHECK( psg->map_dimensions() == 16 + extra );
    }

    auto psg2 = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "sq~uq W+ W+ dq cq~",
        constants
    );

    CHECK( psg2->channels()       == 43 );
    CHECK( psg2->dimensions()     == 10 );
    CHECK( psg2->map_dimensions() == 24 );

    // pp -> 2 jets + two pairs of leptons
    auto psg3 = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "sq~uq ne el~nm mu~dq cq~",
        constants
    );

    CHECK( psg3->channels()       == 93 );
    CHECK( psg3->dimensions()     == 16 );
    CHECK( psg3->map_dimensions() == 32 );

    // pp -> 3 jets + two pairs of leptons
    auto psg4 = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "sq~uq ne el~nm mu~dq cq~gl ",
        constants
    );

    CHECK( psg4->channels()       == 452 );
    CHECK( psg4->dimensions()     == 19 );
    CHECK( psg4->map_dimensions() == 36 );

    // pp -> 2jets

    // NOTE: generator does not generate three-gluon vertices

    // four-gluon vertex
    auto psg5a = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "gl gl gl gl ",
        constants
    );

    CHECK( psg5a->channels()       == 1 );
    CHECK( psg5a->dimensions()     == 4 );
    CHECK( psg5a->map_dimensions() == 16 );

    // t- and u-channel
    auto psg5b = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "uq~uq gl gl ",
        constants
    );

    CHECK( psg5b->channels()       == 2 );
    CHECK( psg5b->dimensions()     == 4 );
    CHECK( psg5b->map_dimensions() == 16 );

    // s-channels with gluon/photon and Z, and t-channel with W
    auto psg5c = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "uq~uq dq dq~",
        constants
    );

    CHECK( psg5c->channels()       == 3 );
    CHECK( psg5c->dimensions()     == 4 );
    CHECK( psg5c->map_dimensions() == 16 );

    // s- and t-channels with gluon/photon and Z
    auto psg5d = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "uq~uq uq uq~",
        constants
    );

    CHECK( psg5d->channels()       == 4 );
    CHECK( psg5d->dimensions()     == 4 );
    CHECK( psg5d->map_dimensions() == 16 );

    auto psg5 = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        std::vector<std::string>{ "gl gl gl gl ", "uq~uq gl gl ",
            "uq~uq dq dq~", "uq~uq uq uq~" },
        constants
    );

    CHECK( psg5->channels()       == 7 );
    CHECK( psg5->dimensions()     == 4 );
    CHECK( psg5->map_dimensions() == 16 );

    // pp -> W+ Z + 2 jets
    auto psg6a = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "uq uq el~ne mu~mu dq uq ",
        constants
    );

    CHECK( psg6a->channels()       == 420 );
    CHECK( psg6a->dimensions()     == 16 );
    CHECK( psg6a->map_dimensions() == 32 );

    // pp -> W+ Z + 2 jets + photon
    auto psg6b = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "uq uq el~ne mu~mu dq uq ga ",
        constants
    );

    CHECK( psg6b->channels()       == 4302 );
    CHECK( psg6b->dimensions()     == 19 );
    CHECK( psg6b->map_dimensions() == 36 );

    // pp -> W+W+W-
    auto psg7 = hep::make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        "dq~uq el ne~mu~nm ta~nt ",
        constants
    );

    CHECK( psg7->channels()       == 85 );
    CHECK( psg7->dimensions()     == 16 );
    CHECK( psg7->map_dimensions() == 32 );
}
