#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"
#include "hep/mc/projector.hpp"

#include "hep/ps/born_integrand.hpp"
#include "hep/ps/final_state.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/lusifer_ps_channels.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/psp.hpp"
#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scales.hpp"
#include "hep/ps/trivial_cutter.hpp"
#include "hep/ps/trivial_recombiner.hpp"

#include "test_phase_space_generator.hpp"

#include "nonstd/span.hpp"

#include "catch2/catch.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <vector>

using T = HEP_TYPE_T;

TEST_CASE("", "[]")
{
    T const energy = T(13000.0);
    T const min_energy = T(20.0);

    T const mw_on = T(80.379);
    T const gw_on = T(2.085);
    T const mz_on = T(91.1876);
    T const gz_on = T(2.4952);

	T const mass_higgs  = T(125.0);
	T const mass_top    = T(173.0);
	T const mass_w      = mw_on / sqrt(T(1.0) + gw_on * gw_on / (mw_on * mw_on));
	T const mass_z      = mz_on / sqrt(T(1.0) + gz_on * gz_on / (mz_on * mz_on));
	T const width_higgs = T(4.088e-3);
	T const width_top   = T();
	T const width_w     = gw_on / sqrt(T(1.0) + gw_on * gw_on / (mw_on * mw_on));
	T const width_z     = gz_on / sqrt(T(1.0) + gz_on * gz_on / (mz_on * mz_on));

    auto const constants = hep::lusifer_constants<T>{mass_higgs, width_higgs, mass_top, width_top,
        mass_w, width_w, mass_z, width_z};

    auto channels = hep::lusifer_ps_channels("uq dq~el ne~mu~nm ta~nt ", constants);

    // channel #60 is the first Higgs channel
    REQUIRE( channels.at(60).invariants().at(0).mom_id() ==  11 );
    REQUIRE( channels.at(60).invariants().at(1).mom_id() ==  47 );
    REQUIRE( channels.at(60).invariants().at(2).mom_id() ==  59 );
    REQUIRE( channels.at(60).invariants().at(3).mom_id() == 191 );
    REQUIRE( channels.at(60).invariants().at(0).pdg_id() ==  24 );
    REQUIRE( channels.at(60).invariants().at(1).pdg_id() ==  24 );
    REQUIRE( channels.at(60).invariants().at(2).pdg_id() ==  25 );
    REQUIRE( channels.at(60).invariants().at(3).pdg_id() ==  24 );

    // channel #63 is the second Higgs channel
    REQUIRE( channels.at(63).invariants().at(0).mom_id() ==  11 );
    REQUIRE( channels.at(63).invariants().at(1).mom_id() ==  47 );
    REQUIRE( channels.at(63).invariants().at(2).mom_id() == 191 );
    REQUIRE( channels.at(63).invariants().at(3).mom_id() == 203 );
    REQUIRE( channels.at(63).invariants().at(0).pdg_id() ==  24 );
    REQUIRE( channels.at(63).invariants().at(1).pdg_id() ==  24 );
    REQUIRE( channels.at(63).invariants().at(2).pdg_id() ==  24 );
    REQUIRE( channels.at(63).invariants().at(3).pdg_id() ==  25 );

    channels.reserve(channels.size() + 4);

    auto channel = channels.at(60);
    std::swap(channel.invariants().at(0), channel.invariants().at(2));
    // 1. Higgs, 2. W- of Higgs decay, 3. W+ of Higgs decay, 4. second W+
    channels.push_back(channel);
    std::swap(channel.invariants().at(1), channel.invariants().at(2));
    // 1. Higgs, 2. W+ of Higgs decay, 3. W- of Higgs decay, 4. second W+
    channels.push_back(channel);
    channel = channels.at(63);
    std::swap(channel.invariants().at(0), channel.invariants().at(1));
    std::swap(channel.invariants().at(1), channel.invariants().at(3));
    //
    channels.push_back(channel);
    std::swap(channel.invariants().at(2), channel.invariants().at(3));
    //
    channels.push_back(channel);

	auto generator = hep::make_lusifer_phase_space_generator(min_energy, energy, channels,
        constants);

    std::vector<T> random_numbers(1000);
    std::vector<T> momenta(4 * 8);
    std::vector<T> densities(89);

    generator->generate(random_numbers, momenta, 60);
    generator->densities(densities);
}
