#include "hep/ps/nfpc.hpp"

#include "catch2/catch.hpp"

using T = HEP_TYPE_T;

TEST_CASE("non-factorizable photonic corrections", "")
{
    T const mw = T(80.385);
    T const gw = T(2.085);

    hep::nfpc_info<T> info;

    info.masses = { mw, mw, mw };
    info.widths = { gw, gw, gw };
    info.non_resonant_particles = { 0, 1 };
    info.decays = { { 2, 3 }, { 4, 5 }, { 6, 7 } };
    info.charges = {
        T( 2.0) / T(3.0), // u
        T(-1.0) / T(3.0), // d
        T(-1.0),          // e
        T(),              // nue
        T(-1.0),          // mu
        T(),              // numu
        T(-1.0),          // tau
        T()               // nutau
    };
    info.signs = {
        T( 1.0), // u
        T(-1.0), // anti-d
        T(-1.0), // e-
        T( 1.0), // anti-nue
        T( 1.0), // mu+
        T(-1.0), // numu
        T( 1.0), // anti-tau
        T(-1.0)  // nutau
    };
    info.photon_mass = mw;

    std::vector<T> const on_shell_momenta = {
        T(1.4036371175956297e+03), T( 0.0000000000000000e+00),
            T( 0.0000000000000000e+00), T( 1.4036371175956297e+03),
        T(1.4036371175956297e+03), T( 0.0000000000000000e+00),
            T( 0.0000000000000000e+00), T(-1.4036371175956297e+03),
        T(4.2861714611646096e+02), T(-3.8145978730431926e+02),
            T(-1.5495928638583607e+02), T( 1.1911636402088367e+02),
        T(3.3584514756097786e+02), T(-3.0414350953916244e+02),
            T(-6.0167589779927802e+01), T( 1.2910673834285620e+02),
        T(3.2941639708202649e+02), T( 3.0475546906530633e+02),
            T( 1.1463344223223046e+00), T(-1.2505179990180332e+02),
        T(1.0540601491900031e+03), T( 9.9887180973483021e+02),
            T( 1.1736949522680268e+02), T(-3.1547156356432833e+02),
        T(1.0197013150371465e+02), T(-8.3696766425332754e+01),
            T( 3.6105764572796822e+01), T( 4.5707031991170922e+01),
        T(5.5736526373807567e+02), T(-5.3432721553132205e+02),
            T( 6.0505281943842107e+01), T( 1.4659322911122092e+02)
    };

    std::vector<T> const virtualities = {
        T(6356.7956255959798), T(6425.0364949135401), T(6311.6500857137944)
    };

    T const result = T(0.0075553042564088764) / acos(T(-1.0)) *
        hep::nfpc(on_shell_momenta, virtualities, info);

    CHECK_THAT( result, Catch::WithinULP(T(1.50109542892689555e-03), 4) );
}
