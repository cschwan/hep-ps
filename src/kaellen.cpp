#include "hep/ps/kaellen.hpp"

#include <cmath>
#include <utility>

namespace
{

// `one` and `two` must not be zero
template <typename T>
bool same_sign_and_relative_equality(T one, T two)
{
	using std::fabs;
	using std::fmax;
	using std::signbit;

	// TODO: what's the best value for `threshold`?
	T const threshold = T(1e-7);

	return (signbit(one) == signbit(two)) &&
		(fabs(one - two) / fmax(fabs(one), fabs(two)) < threshold);
}

}

namespace hep
{

template <typename T>
T kaellen(T x, T y, T z)
{
	using std::copysign;
	using std::fabs;
	using std::swap;

	// sort (x,y,z) in descending order
	if (x < y) { swap(x, y); }
	if (y < z) { swap(y, z); }
	if (x < y) { swap(x, y); }

	if (y == T{})
	{
		if (x == T{})
		{
			return z * z;
		}
		else if (z == T{})
		{
			return x * x;
		}
	}
	else if (same_sign_and_relative_equality(x, y))
	{
		T const sqrt_x = sqrt(fabs(x));
		T const sqrt_y = sqrt(fabs(y));
		T const zsign = copysign(z, y);

		return (zsign - (sqrt_x + sqrt_y) * (sqrt_x + sqrt_y)) *
			(zsign - (sqrt_x - sqrt_y) * (sqrt_x - sqrt_y));
	}
	else if (same_sign_and_relative_equality(y, z))
	{
		T const sqrt_y = sqrt(fabs(y));
		T const sqrt_z = sqrt(fabs(z));
		T const xsign = copysign(x, y);

		return (xsign - (sqrt_y + sqrt_z) * (sqrt_y + sqrt_z)) *
			(xsign - (sqrt_y - sqrt_z) * (sqrt_y - sqrt_z));
	}

	return (x - y - z) * (x - y - z) - T(4.0) * y * z;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template double kaellen(double, double, double);

}
