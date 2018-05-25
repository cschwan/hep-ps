#include "hep/ps/tensor_decompose.hpp"

#include <cmath>

namespace hep
{

template <typename T>
void tensor_decompose(
	std::array<T, 4> const& momentum,
	T denominator,
	std::array<std::array<T, 4>, 4>& results,
	std::array<T, 4>& factors
) {
	using std::copysign;
	using std::fabs;
	using std::sqrt;

	T const p0 = momentum[0];
	T const p1 = momentum[1];
	T const p2 = momentum[2];
	T const p3 = momentum[3];

	T const plen = sqrt(p1 * p1 + p2 * p2 + p3 * p3);
	T const psqr = p0 * p0 - (p1 * p1 + p2 * p2 + p3 * p3);
	T const cos_theta = p3 / plen;
	T const sin_theta = sqrt((T(1.0) - cos_theta) * (T(1.0) + cos_theta));
	T const factor = (T(1.0) - denominator) / denominator / psqr;
	T const x = sqrt(fabs(factor));
	T const sqrt_psqr = sqrt(fabs(psqr));

	T cos_phi;
	T sin_phi;

	if (sin_theta == T())
	{
		cos_phi = T(1.0);
		sin_phi = T();
	}
	else
	{
		cos_phi = p1 / plen / sin_theta;
		sin_phi = p2 / plen / sin_theta;
	}

	results[0] = { T(), cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta };
	results[1] = { T(), sin_phi, -cos_phi, T() };
	results[2] = {
		plen / sqrt_psqr,
		p0 * sin_theta * cos_phi / sqrt_psqr,
		p0 * sin_theta * sin_phi / sqrt_psqr,
		p0 * cos_theta / sqrt_psqr
	};
	results[3] = { x * p0, x * p1, x * p2, x * p3 };

	factors = {
		T(1.0),
		T(1.0),
		copysign(T(1.0), psqr),
		copysign(T(1.0), factor)
	};
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template void tensor_decompose(
	std::array<double, 4> const&,
	double,
	std::array<std::array<double, 4>, 4>&,
	std::array<double, 4>&
);

}
