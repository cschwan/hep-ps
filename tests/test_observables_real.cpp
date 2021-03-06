#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"

#include "hep/ps/initial_state.hpp"
#include "hep/ps/observables_real.hpp"
#include "hep/ps/rambo_phase_space_generator.hpp"
#include "hep/ps/trivial_distributions.hpp"

#include "test_structures.hpp"

#include <catch.hpp>

#include <cstddef>
#include <functional>
#include <vector>

using T = HEP_TYPE_T;

void test_observables_real(
    hep::initial_state_set set,
    std::size_t count,
    std::size_t emitter,
    std::size_t unresolved,
    std::size_t spectator,
    bool inclusive
) {
    auto generator = hep::make_rambo_phase_space_generator<T>(
        T(1.0),
        T(100.0),
        count + 1
    );

    auto integrand = hep::make_observables_real<T>(
        test_matrix_elements<T>(set, count, emitter, unresolved, spectator),
        test_subtraction<T>(emitter, unresolved, spectator),
        test_cuts<T>(count, inclusive),
        test_recombiner<T>(count),
        test_pdf<T>(),
        test_scale_setter<T>(),
        hep::trivial_distributions<T>(),
        set,
        T(1.0),
        T()
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

    if (inclusive)
    {
        CHECK_THAT( result.front().value() ,
            Catch::WithinULP(T(6.3165579739733147e-11), 5) );
    }
    else
    {
        CHECK_THAT( result.front().value() ,
            Catch::WithinULP(T(-6.3165579739733147e-11), 5) );
    }
}

TEST_CASE("observables_real", "[observables_real]")
{
    hep::initial_state_set set{hep::initial_state::q43_cu};

    test_observables_real(set, 4, 2, 3, 0, true);
    test_observables_real(set, 4, 2, 3, 0, false);
}
