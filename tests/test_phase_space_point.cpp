#include "hep/ps/phase_space_point.hpp"
#include "hep/ps/rambo_phase_space_generator.hpp"

#include <catch.hpp>

#include <cmath>
#include <vector>

using T = HEP_TYPE_T;

TEST_CASE("check phase_space_point functions", "[phase_space_point]")
{
    using std::acos;
    using std::atanh;

    T const pi = acos(T(-1.0));
    T const va = atanh(T(1.0) / sqrt(T(3.0)));
    T const ot = T(1.0) / T(3.0);

    std::vector<T> point1 = {
        sqrt(T(3.0)), T( 1.0), T( 1.0), T( 1.0),
        sqrt(T(3.0)), T(-1.0), T( 1.0), T( 1.0),
        sqrt(T(3.0)), T( 1.0), T(-1.0), T( 1.0),
        sqrt(T(3.0)), T( 1.0), T( 1.0), T(-1.0),
        sqrt(T(3.0)), T(-1.0), T(-1.0), T( 1.0),
        sqrt(T(3.0)), T(-1.0), T( 1.0), T(-1.0),
        sqrt(T(3.0)), T( 1.0), T(-1.0), T(-1.0),
        sqrt(T(3.0)), T(-1.0), T(-1.0), T(-1.0),
    };

    hep::phase_space_point<T> ps1{point1};

    CHECK_THAT( ps1.pt2(0) , Catch::WithinULP(T(2.0), 0) );
    CHECK_THAT( ps1.pt2(1) , Catch::WithinULP(T(2.0), 0) );
    CHECK_THAT( ps1.pt2(2) , Catch::WithinULP(T(2.0), 0) );
    CHECK_THAT( ps1.pt2(3) , Catch::WithinULP(T(2.0), 0) );
    CHECK_THAT( ps1.pt2(4) , Catch::WithinULP(T(2.0), 0) );
    CHECK_THAT( ps1.pt2(5) , Catch::WithinULP(T(2.0), 0) );
    CHECK_THAT( ps1.pt2(6) , Catch::WithinULP(T(2.0), 0) );
    CHECK_THAT( ps1.pt2(7) , Catch::WithinULP(T(2.0), 0) );

    CHECK_THAT( ps1.phi(0) , Catch::WithinULP(T( 0.25) * pi, 0) );
    CHECK_THAT( ps1.phi(1) , Catch::WithinULP(T( 0.75) * pi, 0) );
    CHECK_THAT( ps1.phi(2) , Catch::WithinULP(T(-0.25) * pi, 0) );
    CHECK_THAT( ps1.phi(3) , Catch::WithinULP(T( 0.25) * pi, 0) );
    CHECK_THAT( ps1.phi(4) , Catch::WithinULP(T(-0.75) * pi, 0) );
    CHECK_THAT( ps1.phi(5) , Catch::WithinULP(T( 0.75) * pi, 0) );
    CHECK_THAT( ps1.phi(6) , Catch::WithinULP(T(-0.25) * pi, 0) );
    CHECK_THAT( ps1.phi(7) , Catch::WithinULP(T(-0.75) * pi, 0) );

    CHECK_THAT( ps1.rap_neg(0) , Catch::WithinULP(T(-1.0) * va, 0) );
    CHECK_THAT( ps1.rap_neg(1) , Catch::WithinULP(T(-1.0) * va, 0) );
    CHECK_THAT( ps1.rap_neg(2) , Catch::WithinULP(T(-1.0) * va, 0) );
    CHECK_THAT( ps1.rap_neg(3) , Catch::WithinULP(T( 1.0) * va, 0) );
    CHECK_THAT( ps1.rap_neg(4) , Catch::WithinULP(T(-1.0) * va, 0) );
    CHECK_THAT( ps1.rap_neg(5) , Catch::WithinULP(T( 1.0) * va, 0) );
    CHECK_THAT( ps1.rap_neg(6) , Catch::WithinULP(T( 1.0) * va, 0) );
    CHECK_THAT( ps1.rap_neg(7) , Catch::WithinULP(T( 1.0) * va, 0) );

    CHECK_THAT( ps1.rap_pos(0) , Catch::WithinULP(T( 1.0) * va, 0) );
    CHECK_THAT( ps1.rap_pos(1) , Catch::WithinULP(T( 1.0) * va, 0) );
    CHECK_THAT( ps1.rap_pos(2) , Catch::WithinULP(T( 1.0) * va, 0) );
    CHECK_THAT( ps1.rap_pos(3) , Catch::WithinULP(T(-1.0) * va, 0) );
    CHECK_THAT( ps1.rap_pos(4) , Catch::WithinULP(T( 1.0) * va, 0) );
    CHECK_THAT( ps1.rap_pos(5) , Catch::WithinULP(T(-1.0) * va, 0) );
    CHECK_THAT( ps1.rap_pos(6) , Catch::WithinULP(T(-1.0) * va, 0) );
    CHECK_THAT( ps1.rap_pos(7) , Catch::WithinULP(T(-1.0) * va, 0) );

    CHECK_THAT( ps1.prap_neg(0) , Catch::WithinULP(ps1.rap_neg(0), 1) );
    CHECK_THAT( ps1.prap_neg(1) , Catch::WithinULP(ps1.rap_neg(1), 1) );
    CHECK_THAT( ps1.prap_neg(2) , Catch::WithinULP(ps1.rap_neg(2), 1) );
    CHECK_THAT( ps1.prap_neg(3) , Catch::WithinULP(ps1.rap_neg(3), 1) );
    CHECK_THAT( ps1.prap_neg(4) , Catch::WithinULP(ps1.rap_neg(4), 1) );
    CHECK_THAT( ps1.prap_neg(5) , Catch::WithinULP(ps1.rap_neg(5), 1) );
    CHECK_THAT( ps1.prap_neg(6) , Catch::WithinULP(ps1.rap_neg(6), 1) );
    CHECK_THAT( ps1.prap_neg(7) , Catch::WithinULP(ps1.rap_neg(7), 1) );

    CHECK_THAT( ps1.prap_pos(0) , Catch::WithinULP(ps1.rap_pos(0), 1) );
    CHECK_THAT( ps1.prap_pos(1) , Catch::WithinULP(ps1.rap_pos(1), 1) );
    CHECK_THAT( ps1.prap_pos(2) , Catch::WithinULP(ps1.rap_pos(2), 1) );
    CHECK_THAT( ps1.prap_pos(3) , Catch::WithinULP(ps1.rap_pos(3), 1) );
    CHECK_THAT( ps1.prap_pos(4) , Catch::WithinULP(ps1.rap_pos(4), 1) );
    CHECK_THAT( ps1.prap_pos(5) , Catch::WithinULP(ps1.rap_pos(5), 1) );
    CHECK_THAT( ps1.prap_pos(6) , Catch::WithinULP(ps1.rap_pos(6), 1) );
    CHECK_THAT( ps1.prap_pos(7) , Catch::WithinULP(ps1.rap_pos(7), 1) );

    CHECK_THAT( ps1.cos_angle_neg(0, 1) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(0, 2) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(0, 3) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(0, 4) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(0, 5) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(0, 6) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(0, 7) , Catch::WithinULP(T(-1.0), 5) );
    CHECK_THAT( ps1.cos_angle_neg(1, 2) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(1, 3) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(1, 4) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(1, 5) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(1, 6) , Catch::WithinULP(T(-1.0), 0) );
    CHECK_THAT( ps1.cos_angle_neg(1, 7) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(2, 3) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(2, 4) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(2, 5) , Catch::WithinULP(T(-1.0), 0) );
    CHECK_THAT( ps1.cos_angle_neg(2, 6) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(2, 7) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(3, 4) , Catch::WithinULP(T(-1.0), 0) );
    CHECK_THAT( ps1.cos_angle_neg(3, 5) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(3, 6) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(3, 7) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(4, 5) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(4, 6) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(4, 7) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(5, 6) , Catch::WithinULP(T(-1.0) * ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(5, 7) , Catch::WithinULP(ot, 5) );
    CHECK_THAT( ps1.cos_angle_neg(6, 7) , Catch::WithinULP(ot, 5) );

    CHECK_THAT( ps1.cos_angle_pos(0, 1) ,
        Catch::WithinULP(ps1.cos_angle_neg(0, 1), 0) );
    CHECK_THAT( ps1.cos_angle_pos(0, 2) ,
        Catch::WithinULP(ps1.cos_angle_neg(0, 2), 0) );
    CHECK_THAT( ps1.cos_angle_pos(0, 3) ,
        Catch::WithinULP(ps1.cos_angle_neg(0, 3), 0) );
    CHECK_THAT( ps1.cos_angle_pos(0, 4) ,
        Catch::WithinULP(ps1.cos_angle_neg(0, 4), 0) );
    CHECK_THAT( ps1.cos_angle_pos(0, 5) ,
        Catch::WithinULP(ps1.cos_angle_neg(0, 5), 0) );
    CHECK_THAT( ps1.cos_angle_pos(0, 6) ,
        Catch::WithinULP(ps1.cos_angle_neg(0, 6), 0) );
    CHECK_THAT( ps1.cos_angle_pos(0, 7) ,
        Catch::WithinULP(ps1.cos_angle_neg(0, 7), 0) );
    CHECK_THAT( ps1.cos_angle_pos(1, 2) ,
        Catch::WithinULP(ps1.cos_angle_neg(1, 2), 0) );
    CHECK_THAT( ps1.cos_angle_pos(1, 3) ,
        Catch::WithinULP(ps1.cos_angle_neg(1, 3), 0) );
    CHECK_THAT( ps1.cos_angle_pos(1, 4) ,
        Catch::WithinULP(ps1.cos_angle_neg(1, 4), 0) );
    CHECK_THAT( ps1.cos_angle_pos(1, 5) ,
        Catch::WithinULP(ps1.cos_angle_neg(1, 5), 0) );
    CHECK_THAT( ps1.cos_angle_pos(1, 6) ,
        Catch::WithinULP(ps1.cos_angle_neg(1, 6), 0) );
    CHECK_THAT( ps1.cos_angle_pos(1, 7) ,
        Catch::WithinULP(ps1.cos_angle_neg(1, 7), 0) );
    CHECK_THAT( ps1.cos_angle_pos(2, 3) ,
        Catch::WithinULP(ps1.cos_angle_neg(2, 3), 0) );
    CHECK_THAT( ps1.cos_angle_pos(2, 4) ,
        Catch::WithinULP(ps1.cos_angle_neg(2, 4), 0) );
    CHECK_THAT( ps1.cos_angle_pos(2, 5) ,
        Catch::WithinULP(ps1.cos_angle_neg(2, 5), 0) );
    CHECK_THAT( ps1.cos_angle_pos(2, 6) ,
        Catch::WithinULP(ps1.cos_angle_neg(2, 6), 0) );
    CHECK_THAT( ps1.cos_angle_pos(2, 7) ,
        Catch::WithinULP(ps1.cos_angle_neg(2, 7), 0) );
    CHECK_THAT( ps1.cos_angle_pos(3, 4) ,
        Catch::WithinULP(ps1.cos_angle_neg(3, 4), 0) );
    CHECK_THAT( ps1.cos_angle_pos(3, 5) ,
        Catch::WithinULP(ps1.cos_angle_neg(3, 5), 0) );
    CHECK_THAT( ps1.cos_angle_pos(3, 6) ,
        Catch::WithinULP(ps1.cos_angle_neg(3, 6), 0) );
    CHECK_THAT( ps1.cos_angle_pos(3, 7) ,
        Catch::WithinULP(ps1.cos_angle_neg(3, 7), 0) );
    CHECK_THAT( ps1.cos_angle_pos(4, 5) ,
        Catch::WithinULP(ps1.cos_angle_neg(4, 5), 0) );
    CHECK_THAT( ps1.cos_angle_pos(4, 6) ,
        Catch::WithinULP(ps1.cos_angle_neg(4, 6), 0) );
    CHECK_THAT( ps1.cos_angle_pos(4, 7) ,
        Catch::WithinULP(ps1.cos_angle_neg(4, 7), 0) );
    CHECK_THAT( ps1.cos_angle_pos(5, 6) ,
        Catch::WithinULP(ps1.cos_angle_neg(5, 6), 0) );
    CHECK_THAT( ps1.cos_angle_pos(5, 7) ,
        Catch::WithinULP(ps1.cos_angle_neg(5, 7), 0) );
    CHECK_THAT( ps1.cos_angle_pos(6, 7) ,
        Catch::WithinULP(ps1.cos_angle_neg(6, 7), 0) );

    // TODO: write checks for
    // - abs_phi_diff
    // - dist2
    // - m2 (2x)
    // - mt (2x)
    // - pt2
    // - rap_diff
}
