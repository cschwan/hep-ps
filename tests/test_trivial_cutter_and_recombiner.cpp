#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"

#include "hep/ps/born_integrand.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/rambo_phase_space_generator.hpp"
#include "hep/ps/trivial_cutter.hpp"
#include "hep/ps/trivial_distributions.hpp"
#include "hep/ps/trivial_recombiner.hpp"

#include "test_structures.hpp"

#include <catch.hpp>

#include <cstddef>
#include <functional>
#include <vector>

using T = HEP_TYPE_T;

void test_trivial_cutter_and_recombiner(
    hep::initial_state_set set,
    std::size_t count
) {
    auto generator = hep::make_rambo_phase_space_generator<T>(
        T(1.0),
        T(100.0),
        count
    );

    auto integrand = hep::make_born_integrand<T>(
        test_matrix_elements<T>(set, count),
        hep::trivial_cutter<T>(),
        hep::trivial_recombiner<T>(),
        test_pdf<T>(),
        test_scale_setter<T>{false},
        hep::trivial_distributions<T>(),
        set,
        T(1.0)
    );

    auto const result = hep::multi_channel(
        hep::make_multi_channel_integrand<T>(
            std::ref(*integrand),
            generator->dimensions(),
            std::ref(*generator),
            generator->map_dimensions(),
            generator->channels(),
            test_distribution_parameters<T>()
        ),
        std::vector<std::size_t>{1}
    );

    CHECK_THAT( result.front().value() ,
        Catch::WithinULP(T(0.0015980493564892501), 5) );
}

TEST_CASE("check trivial cutter and recombiner", "[trivial_cutter],[trivial_recombiner]")
{
    hep::initial_state_set set{hep::initial_state::cq_uq};

    test_trivial_cutter_and_recombiner(set, 4);
}
