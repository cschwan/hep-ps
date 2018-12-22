#include "hep/ps/ps_functions.hpp"

#include <cmath>

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

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template void boost(double, std::array<double, 4>&, std::array<double, 4> const&, bool);

template void rotate(std::array<double, 4>&, double, double);

}
