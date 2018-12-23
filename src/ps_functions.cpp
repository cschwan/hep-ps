#include "hep/ps/ps_functions.hpp"

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
void boost(T m, std::array<T, 4>& p, std::array<T, 4> const& q, bool inverse)
{
    T const sign = inverse ? T(1.0) : T(-1.0);
    T const bx = sign * q[1] / m;
    T const by = sign * q[2] / m;
    T const bz = sign * q[3] / m;
    T const gamma = q[0] / m;
    T const a = T(1.0) / (T(1.0) + gamma);
    T const p0 = p[0];
    T const px = p[1];
    T const py = p[2];
    T const pz = p[3];
    T const bp = bx * px + by * py + bz * pz;

    p[0] = gamma * p0 + bp;
    p[1] = px + bx * p0 + a * bp * bx;
    p[2] = py + by * p0 + a * bp * by;
    p[3] = pz + bz * p0 + a * bp * bz;
}

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

    if (y == T())
    {
        if (x == T())
        {
            return z * z;
        }
        else if (z == T())
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

template <typename T>
T sqrt_kaellen(T x, T y, T z)
{
    using std::fabs;
    using std::sqrt;

    if (x == T())
    {
        return fabs(y + z);
    }
    else if (y == T())
    {
        return fabs(x - z);
    }
    else if (z == T())
    {
        return fabs(x - y);
    }

    T const diff = fabs(x - y - z);
    T const arg = T(4.0) * y * z / (diff * diff);
    T const result = diff * sqrt(fabs(T(1.0) - arg));

    return result;
}

template <typename T>
void rotate(std::array<T, 4>& p, T phi, T cos_theta)
{
    using std::cos;
    using std::fabs;
    using std::sin;
    using std::sqrt;

    T const sin_phi = sin(phi);
    T const cos_phi = cos(phi);
    T const sin_theta = sqrt(fabs((T(1.0) - cos_theta) * (T(1.0) + cos_theta)));

    T const px = p[1];
    T const py = p[2];
    T const pz = p[3];

    p[1] =  px * cos_theta * cos_phi + py * sin_phi + pz * sin_theta * cos_phi;
    p[2] = -px * cos_theta * sin_phi + py * cos_phi - pz * sin_theta * sin_phi;
    p[3] = -px * sin_theta                          + pz * cos_theta;
}

template <typename T>
space_like_invariant_bounds<T> calculate_space_like_invariant_bounds(T s, T s1, T s2, T t1, T t2)
{
    using std::fabs;
    using std::sqrt;

    T const threshold = T(1e-5);

    std::size_t const non_zero_invariants =
        ((s1 == T()) ? 0 : 1) | ((s2 == T()) ? 0 : 2) |
        ((t1 == T()) ? 0 : 4) | ((t2 == T()) ? 0 : 8);

    T const lambdas = sqrt_kaellen(s, s1, s2);
    T const lambdat = sqrt_kaellen(s, t1, t2);

    T tmin = T();
    T tmax = T();

    T const x = s - s1 - s2;
    T const y = s - t1 - t2;
    T const e1 = s1 * s2 / (x * x);
    T const e2 = t1 * t2 / (y * y);

    switch (non_zero_invariants)
    {
    case 1: // s2 = t1 = t2 = 0
    case 2: // s1 = t1 = t2 = 0
        tmin = -lambdas;
        // tmax is zero
        break;

    case 3: // t1 = t2 = 0
        tmin = T(-0.5) * x * (T(1.0) + sqrt(fabs(T(1.0) - T(4.0) * e1)));
        tmax = T(-0.5) * x * (T(1.0) - sqrt(fabs(T(1.0) - T(4.0) * e1)));
        // TODO: separate code path for small e1
        break;

    case 4: // s1 = s2 = t2 = 0
    case 8: // s1 = s2 = t1 = 0
        tmin = -lambdat;
        // tmax is zero
        break;

    case 5: // s2 = t2 = 0
        tmin = -s + s1 + t1 - s1 * t1 / s;
        // tmax is zero
        break;

    case 6: // s1 = t2 = 0
        tmin = -s + s2 + t1;
        tmax = s2 * t1 / s;
        break;

    case 9: // s2 = t1 = 0
        tmin = -s + s1 + t2;
        tmax = s1 * t2 / s;
        break;

    case 10: // s1 = t1 = 0
        tmin = -s + s2 + t2 - s2 * t2 / s;
        // tmax is zero
        break;

    case 12: // s1 = s2 = 0
        tmin = T(-0.5) * y * (T(1.0) + sqrt(fabs(T(1.0) - T(4.0) * e2)));
        tmax = T(-0.5) * y * (T(1.0) - sqrt(fabs(T(1.0) - T(4.0) * e2)));
        // TODO: separate code path for small e2
        break;

    default:
        T tmp = (s + s1 - s2) * (s + t1 - t2);
        tmin = s1 + t1 - T(0.5) * (tmp + lambdas * lambdat) / s;
        tmax = s1 + t1 - T(0.5) * (tmp - lambdas * lambdat) / s;

        if ((e1 < threshold) && (e2 < threshold))
        {
            // TODO: add order six terms?
            tmax = ((t2 * s1 + t1 * s2) - (x * y * (e1 + e2 + (e1 - e2) * (e1 - e2)))) / s;
        }
    }

    return { tmin, tmax, lambdas, lambdat };
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template void boost(double, std::array<double, 4>&, std::array<double, 4> const&, bool);

template double kaellen(double, double, double);

template double sqrt_kaellen(double, double, double);

template void rotate(std::array<double, 4>&, double, double);

template space_like_invariant_bounds<double>
calculate_space_like_invariant_bounds(double, double, double, double, double);

}
