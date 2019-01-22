#include "hep/ps/psp.hpp"

#include <cmath>
#include <iterator>

namespace hep
{

template <typename T>
psp<T>::psp(
    std::vector<T>& cms_psp,
    std::vector<recombined_state>& states,
    T rap_shift,
    psp_type type
)
    : p_{cms_psp}
    , states_{states}
    , rap_shift_{rap_shift}
    , sign_{(type == psp_type::pos_rap) ? T(1.0) : T(-1.0)}
{
}

template <typename T>
T psp<T>::abs_phi_diff(std::size_t i, std::size_t j) const
{
    using std::acos;
    using std::fabs;

    T const result = fabs(phi(i) - phi(j));
    T const pi = acos(T(-1.0));

    return (result > pi) ? (T(2.0) * pi - result) : result;
}

template <typename T>
T psp<T>::cos_angle(std::size_t i, std::size_t j) const
{
    using std::cosh;
    using std::sinh;
    using std::sqrt;

    T const px = p_[4*i+1];
    T const py = p_[4*i+2];
    T const qx = p_[4*j+1];
    T const qy = p_[4*j+2];

    T const x1 = p_[4*i+3] / p_[4*i+0];
    T const x2 = p_[4*j+3] / p_[4*j+0];
    T const cosh_rap = cosh(rap_shift_);
    T const sinh_rap = sinh(rap_shift_);
    T const pt = sqrt(px * px + py * py);
    T const qt = sqrt(qx * qx + qy * qy);

    T const pz = pt * (sinh_rap + sign_ * cosh_rap * x1) / sqrt(T(1.0) - x1 * x1);
    T const qz = qt * (sinh_rap + sign_ * cosh_rap * x2) / sqrt(T(1.0) - x2 * x2);

    T const num = px*qx + py*qy + pz*qz;
    T const den = sqrt((px*px + py*py + pz*pz) * (qx*qx + qy*qy + qz*qz));

    return num / den;
}

template <typename T>
T psp<T>::dist2(std::size_t i, std::size_t j) const
{
    T const rap = rap_diff(i, j);
    T const phi = abs_phi_diff(i, j);

    return rap * rap + phi * phi;
}

template <typename T>
T psp<T>::m2(std::size_t i) const
{
    T const p0 = p_[4*i+0];
    T const p1 = p_[4*i+1];
    T const p2 = p_[4*i+2];
    T const p3 = p_[4*i+3];

    return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

template <typename T>
T psp<T>::m2(std::size_t i, std::size_t j) const
{
    T const p0 = p_[4*i+0] + p_[4*j+0];
    T const p1 = p_[4*i+1] + p_[4*j+1];
    T const p2 = p_[4*i+2] + p_[4*j+2];
    T const p3 = p_[4*i+3] + p_[4*j+3];

    return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

template <typename T>
T psp<T>::m2(std::size_t i, std::size_t j, std::size_t k) const
{
    T const p0 = p_[4*i+0] + p_[4*j+0] + p_[4*k+0];
    T const p1 = p_[4*i+1] + p_[4*j+1] + p_[4*k+1];
    T const p2 = p_[4*i+2] + p_[4*j+2] + p_[4*k+2];
    T const p3 = p_[4*i+3] + p_[4*j+3] + p_[4*k+3];

    return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

template <typename T>
T psp<T>::m2(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const
{
    T const p0 = p_[4*i+0] + p_[4*j+0] + p_[4*k+0] + p_[4*l+0];
    T const p1 = p_[4*i+1] + p_[4*j+1] + p_[4*k+1] + p_[4*l+1];
    T const p2 = p_[4*i+2] + p_[4*j+2] + p_[4*k+2] + p_[4*l+2];
    T const p3 = p_[4*i+3] + p_[4*j+3] + p_[4*k+3] + p_[4*l+3];

    return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

template <typename T>
T psp<T>::m2(std::size_t i, std::size_t j, std::size_t k, std::size_t l, std::size_t m) const
{
    T const p0 = p_[4*i+0] + p_[4*j+0] + p_[4*k+0] + p_[4*l+0] + p_[4*m+0];
    T const p1 = p_[4*i+1] + p_[4*j+1] + p_[4*k+1] + p_[4*l+1] + p_[4*m+1];
    T const p2 = p_[4*i+2] + p_[4*j+2] + p_[4*k+2] + p_[4*l+2] + p_[4*m+2];
    T const p3 = p_[4*i+3] + p_[4*j+3] + p_[4*k+3] + p_[4*l+3] + p_[4*m+3];

    return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

template <typename T>
T psp<T>::m2(std::size_t i, std::size_t j, std::size_t k, std::size_t l, std::size_t m, std::size_t n) const
{
    T const p0 = p_[4*i+0] + p_[4*j+0] + p_[4*k+0] + p_[4*l+0] + p_[4*m+0] + p_[4*n+0];
    T const p1 = p_[4*i+1] + p_[4*j+1] + p_[4*k+1] + p_[4*l+1] + p_[4*m+1] + p_[4*n+1];
    T const p2 = p_[4*i+2] + p_[4*j+2] + p_[4*k+2] + p_[4*l+2] + p_[4*m+2] + p_[4*n+2];
    T const p3 = p_[4*i+3] + p_[4*j+3] + p_[4*k+3] + p_[4*l+3] + p_[4*m+3] + p_[4*n+3];

    return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

template <typename T>
T psp<T>::mt(std::size_t i, std::size_t j) const
{
    using std::sqrt;

    T const pt = sqrt(pt2(i)) + sqrt(pt2(j));
    T const px = p_[4*i+1] + p_[4*j+1];
    T const py = p_[4*i+2] + p_[4*j+2];

    return sqrt(pt * pt - px * px - py * py);
}

template <typename T>
T psp<T>::mt(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const
{
    using std::sqrt;

    T const pt = sqrt(pt2(i))+sqrt(pt2(j))+sqrt(pt2(k))+sqrt(pt2(l));
    T const px = p_[4*i+1] + p_[4*j+1] + p_[4*k+1] + p_[4*l+1];
    T const py = p_[4*i+2] + p_[4*j+2] + p_[4*k+2] + p_[4*l+2];

    return sqrt(pt * pt - px * px - py * py);
}

template <typename T>
T psp<T>::pdist2(std::size_t i, std::size_t j) const
{
    T const rap = prap_diff(i, j);
    T const phi = abs_phi_diff(i, j);

    return rap * rap + phi * phi;
}

template <typename T>
T psp<T>::phi(std::size_t i) const
{
    using std::atan2;

    T const p1 = p_[4*i+1];
    T const p2 = p_[4*i+2];

    return atan2(p2, p1);
}

template <typename T>
T psp<T>::pt2(std::size_t i) const
{
    T const p1 = p_[4*i+1];
    T const p2 = p_[4*i+2];

    return p1 * p1 + p2 * p2;
}

template <typename T>
T psp<T>::pt2(std::size_t i, std::size_t j) const
{
    T const p1 = p_[4*i+1] + p_[4*j+1];
    T const p2 = p_[4*i+2] + p_[4*j+2];

    return p1 * p1 + p2 * p2;
}

template <typename T>
T psp<T>::pt2(std::size_t i, std::size_t j, std::size_t k) const
{
    T const p1 = p_[4*i+1] + p_[4*j+1] + p_[4*k+1];
    T const p2 = p_[4*i+2] + p_[4*j+2] + p_[4*k+2];

    return p1 * p1 + p2 * p2;
}

template <typename T>
T psp<T>::rap_diff(std::size_t i, std::size_t j) const
{
    using std::atanh;

    T const x = p_[4*i+3] / p_[4*i+0];
    T const y = p_[4*j+3] / p_[4*j+0];

    return sign_ * atanh((x - y) / (T(1.0) - x * y));
}

template <typename T>
T psp<T>::prap_diff(std::size_t i, std::size_t j) const
{
    using std::fabs;

    return fabs(prap(i) - prap(j));
}

template <typename T>
T psp<T>::rap(std::size_t i) const
{
    using std::atanh;

    T const p0 = p_[4*i+0];
    T const p3 = p_[4*i+3];

    return rap_shift_ + sign_ * atanh(p3 / p0);
}

template <typename T>
T psp<T>::prap(std::size_t i) const
{
    using std::atanh;
    using std::sinh;
    using std::sqrt;

    T const p0 = p_[4*i+0];
    T const p1 = p_[4*i+1];
    T const p2 = p_[4*i+2];
    T const p3 = p_[4*i+3];

    T const mt = sqrt((p0 - p3) * (p0 + p3));
    T const pt2 = p1 * p1 + p2 * p2;
    T const sinh_rap = sinh(rap_shift_ + sign_ * atanh(p3 / p0));
    T const arg = mt * sinh_rap / sqrt(pt2 + mt * mt * sinh_rap * sinh_rap);

    return atanh(arg);
}

template <typename T>
std::vector<T> const& psp<T>::cms_psp() const
{
    return p_;
}

template <typename T>
T psp<T>::rap_shift() const
{
    return rap_shift_;
}

template <typename T>
psp_type psp<T>::type() const
{
    return (sign_ == T(1.0)) ? psp_type::pos_rap : psp_type::neg_rap;
}

template <typename T>
void psp<T>::remove(std::size_t i)
{
    p_.erase(std::next(p_.begin(), 4 * i), std::next(p_.begin(), 4 * (i + 1)));
    states_.erase(std::next(states_.begin(), i - 2));
}

template <typename T>
recombined_state psp<T>::state(std::size_t i) const
{
    return states_.at(i);
}

template <typename T>
std::vector<recombined_state> const& psp<T>::states() const
{
    return states_;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class psp<double>;

}
