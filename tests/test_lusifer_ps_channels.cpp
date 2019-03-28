#include "hep/ps/lusifer_constants.hpp"
#include "hep/ps/lusifer_ps_channels.hpp"

#include "catch2/catch.hpp"

#include <string>
#include <vector>

using T = double;

hep::lusifer_constants<T> constants(
    T(125.09), T(4.0e-3),
    T(174.2), T(1.41),
    T(80.385), T(2.085),
    T(91.1876), T(2.4952)
);

TEST_CASE("constructors", "[lusifer_phase_space_generator]")
{
    auto channels1 = hep::lusifer_ps_channels("el~el mu mu~", constants);

    // two s-channels with a photon or Z boson
    REQUIRE( channels1.size() == 2 );

    // no invariants
    CHECK( channels1.at(0).invariants().empty() );
    CHECK( channels1.at(1).invariants().empty() );

    // no t-channels
    CHECK( channels1.at(0).tchannels().empty() );
    CHECK( channels1.at(1).tchannels().empty() );

    // one decay for each channel
    CHECK( channels1.at(0).decays().size() == 1 );
    CHECK( channels1.at(1).decays().size() == 1 );

    auto channels2 = hep::lusifer_ps_channels("sq~uq W+ W+ dq cq~", constants);

    CHECK( channels2.size() == 43 );

    // pp -> 2 jets + two pairs of leptons
    auto channels3 = hep::lusifer_ps_channels("sq~uq ne el~nm mu~dq cq~", constants);

    CHECK( channels3.size() == 93 );

    // pp -> 3 jets + two pairs of leptons
    auto channels4 = hep::lusifer_ps_channels("sq~uq ne el~nm mu~dq cq~gl ", constants);

    CHECK( channels4.size() == 452 );

    // pp -> 2jets

    // NOTE: generator does not generate three-gluon vertices

    // four-gluon vertex
    auto channels5a = hep::lusifer_ps_channels("gl gl gl gl ", constants);

    CHECK( channels5a.size() == 1 );

    // t- and u-channel
    auto channels5b = hep::lusifer_ps_channels("uq~uq gl gl ", constants);

    CHECK( channels5b.size() == 2 );

    // s-channels with gluon/photon and Z, and t-channel with W
    auto channels5c = hep::lusifer_ps_channels("uq~uq dq dq~", constants);

    CHECK( channels5c.size() == 3 );

    // s- and t-channels with gluon/photon and Z
    auto channels5d = hep::lusifer_ps_channels("uq~uq uq uq~", constants);

    CHECK( channels5d.size() == 4 );

    // everything together
    auto channels5 = hep::lusifer_ps_channels(
        std::vector<std::string>{ "gl gl gl gl ", "uq~uq gl gl ", "uq~uq dq dq~", "uq~uq uq uq~" },
        constants
    );

    CHECK( channels5.size() == 7 );

    // pp -> W+ Z + 2 jets
    auto channels6a = hep::lusifer_ps_channels("uq uq el~ne mu~mu dq uq ", constants);

    CHECK( channels6a.size() == 420 );

    // pp -> W+ Z + 2 jets + photon
    auto channels6b = hep::lusifer_ps_channels("uq uq el~ne mu~mu dq uq ga ", constants);

    CHECK( channels6b.size() == 4302 );

    // pp -> W+W+W-
    auto channels7 = hep::lusifer_ps_channels("dq~uq el ne~mu~nm ta~nt ", constants);

    CHECK( channels7.size() == 85 );
}
