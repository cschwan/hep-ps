#include "hep/ps/phase_space_point.hpp"
#include "hep/ps/rambo_phase_space_generator.hpp"

#include "catch.hpp"

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

	CHECK( ps1.pt2(0) == T(2.0) );
	CHECK( ps1.pt2(1) == T(2.0) );
	CHECK( ps1.pt2(2) == T(2.0) );
	CHECK( ps1.pt2(3) == T(2.0) );
	CHECK( ps1.pt2(4) == T(2.0) );
	CHECK( ps1.pt2(5) == T(2.0) );
	CHECK( ps1.pt2(6) == T(2.0) );
	CHECK( ps1.pt2(7) == T(2.0) );

	CHECK( ps1.phi(0) == T( 0.25) * pi );
	CHECK( ps1.phi(1) == T( 0.75) * pi );
	CHECK( ps1.phi(2) == T(-0.25) * pi );
	CHECK( ps1.phi(3) == T( 0.25) * pi );
	CHECK( ps1.phi(4) == T(-0.75) * pi );
	CHECK( ps1.phi(5) == T( 0.75) * pi );
	CHECK( ps1.phi(6) == T(-0.25) * pi );
	CHECK( ps1.phi(7) == T(-0.75) * pi );

	CHECK( ps1.rap_neg(0) == T(-1.0) * va );
	CHECK( ps1.rap_neg(1) == T(-1.0) * va );
	CHECK( ps1.rap_neg(2) == T(-1.0) * va );
	CHECK( ps1.rap_neg(3) == T( 1.0) * va );
	CHECK( ps1.rap_neg(4) == T(-1.0) * va );
	CHECK( ps1.rap_neg(5) == T( 1.0) * va );
	CHECK( ps1.rap_neg(6) == T( 1.0) * va );
	CHECK( ps1.rap_neg(7) == T( 1.0) * va );

	CHECK( ps1.rap_pos(0) == T( 1.0) * va );
	CHECK( ps1.rap_pos(1) == T( 1.0) * va );
	CHECK( ps1.rap_pos(2) == T( 1.0) * va );
	CHECK( ps1.rap_pos(3) == T(-1.0) * va );
	CHECK( ps1.rap_pos(4) == T( 1.0) * va );
	CHECK( ps1.rap_pos(5) == T(-1.0) * va );
	CHECK( ps1.rap_pos(6) == T(-1.0) * va );
	CHECK( ps1.rap_pos(7) == T(-1.0) * va );

	CHECK( ps1.prap_neg(0) == Approx(ps1.rap_neg(0)) );
	CHECK( ps1.prap_neg(1) == Approx(ps1.rap_neg(1)) );
	CHECK( ps1.prap_neg(2) == Approx(ps1.rap_neg(2)) );
	CHECK( ps1.prap_neg(3) == Approx(ps1.rap_neg(3)) );
	CHECK( ps1.prap_neg(4) == Approx(ps1.rap_neg(4)) );
	CHECK( ps1.prap_neg(5) == Approx(ps1.rap_neg(5)) );
	CHECK( ps1.prap_neg(6) == Approx(ps1.rap_neg(6)) );
	CHECK( ps1.prap_neg(7) == Approx(ps1.rap_neg(7)) );

	CHECK( ps1.prap_pos(0) == Approx(ps1.rap_pos(0)) );
	CHECK( ps1.prap_pos(1) == Approx(ps1.rap_pos(1)) );
	CHECK( ps1.prap_pos(2) == Approx(ps1.rap_pos(2)) );
	CHECK( ps1.prap_pos(3) == Approx(ps1.rap_pos(3)) );
	CHECK( ps1.prap_pos(4) == Approx(ps1.rap_pos(4)) );
	CHECK( ps1.prap_pos(5) == Approx(ps1.rap_pos(5)) );
	CHECK( ps1.prap_pos(6) == Approx(ps1.rap_pos(6)) );
	CHECK( ps1.prap_pos(7) == Approx(ps1.rap_pos(7)) );

	CHECK( ps1.cos_angle_neg(0, 1) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(0, 2) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(0, 3) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(0, 4) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(0, 5) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(0, 6) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(0, 7) == T(-1.0) );
	CHECK( ps1.cos_angle_neg(1, 2) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(1, 3) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(1, 4) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(1, 5) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(1, 6) == T(-1.0) );
	CHECK( ps1.cos_angle_neg(1, 7) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(2, 3) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(2, 4) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(2, 5) == T(-1.0) );
	CHECK( ps1.cos_angle_neg(2, 6) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(2, 7) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(3, 4) == T(-1.0) );
	CHECK( ps1.cos_angle_neg(3, 5) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(3, 6) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(3, 7) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(4, 5) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(4, 6) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(4, 7) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(5, 6) == Approx(T(-1.0) * ot) );
	CHECK( ps1.cos_angle_neg(5, 7) == Approx(ot) );
	CHECK( ps1.cos_angle_neg(6, 7) == Approx(ot) );

	CHECK( ps1.cos_angle_pos(0, 1) == ps1.cos_angle_neg(0, 1) );
	CHECK( ps1.cos_angle_pos(0, 2) == ps1.cos_angle_neg(0, 2) );
	CHECK( ps1.cos_angle_pos(0, 3) == ps1.cos_angle_neg(0, 3) );
	CHECK( ps1.cos_angle_pos(0, 4) == ps1.cos_angle_neg(0, 4) );
	CHECK( ps1.cos_angle_pos(0, 5) == ps1.cos_angle_neg(0, 5) );
	CHECK( ps1.cos_angle_pos(0, 6) == ps1.cos_angle_neg(0, 6) );
	CHECK( ps1.cos_angle_pos(0, 7) == ps1.cos_angle_neg(0, 7) );
	CHECK( ps1.cos_angle_pos(1, 2) == ps1.cos_angle_neg(1, 2) );
	CHECK( ps1.cos_angle_pos(1, 3) == ps1.cos_angle_neg(1, 3) );
	CHECK( ps1.cos_angle_pos(1, 4) == ps1.cos_angle_neg(1, 4) );
	CHECK( ps1.cos_angle_pos(1, 5) == ps1.cos_angle_neg(1, 5) );
	CHECK( ps1.cos_angle_pos(1, 6) == ps1.cos_angle_neg(1, 6) );
	CHECK( ps1.cos_angle_pos(1, 7) == ps1.cos_angle_neg(1, 7) );
	CHECK( ps1.cos_angle_pos(2, 3) == ps1.cos_angle_neg(2, 3) );
	CHECK( ps1.cos_angle_pos(2, 4) == ps1.cos_angle_neg(2, 4) );
	CHECK( ps1.cos_angle_pos(2, 5) == ps1.cos_angle_neg(2, 5) );
	CHECK( ps1.cos_angle_pos(2, 6) == ps1.cos_angle_neg(2, 6) );
	CHECK( ps1.cos_angle_pos(2, 7) == ps1.cos_angle_neg(2, 7) );
	CHECK( ps1.cos_angle_pos(3, 4) == ps1.cos_angle_neg(3, 4) );
	CHECK( ps1.cos_angle_pos(3, 5) == ps1.cos_angle_neg(3, 5) );
	CHECK( ps1.cos_angle_pos(3, 6) == ps1.cos_angle_neg(3, 6) );
	CHECK( ps1.cos_angle_pos(3, 7) == ps1.cos_angle_neg(3, 7) );
	CHECK( ps1.cos_angle_pos(4, 5) == ps1.cos_angle_neg(4, 5) );
	CHECK( ps1.cos_angle_pos(4, 6) == ps1.cos_angle_neg(4, 6) );
	CHECK( ps1.cos_angle_pos(4, 7) == ps1.cos_angle_neg(4, 7) );
	CHECK( ps1.cos_angle_pos(5, 6) == ps1.cos_angle_neg(5, 6) );
	CHECK( ps1.cos_angle_pos(5, 7) == ps1.cos_angle_neg(5, 7) );
	CHECK( ps1.cos_angle_pos(6, 7) == ps1.cos_angle_neg(6, 7) );

	// TODO: write checks for
	// - abs_phi_diff
	// - dist2
	// - m2 (2x)
	// - mt (2x)
	// - pt2
	// - rap_diff
}
