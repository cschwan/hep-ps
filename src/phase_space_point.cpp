#include "hep/ps/phase_space_point.hpp"

#include <cassert>
#include <cmath>

namespace hep
{

template <typename T>
phase_space_point<T>::phase_space_point(
	std::vector<T> const& phase_space,
	T rapidity_shift
)
	: p_{phase_space.data()}
	, rapidity_shift_{rapidity_shift}
{
}

template <typename T>
T phase_space_point<T>::abs_phi_diff(std::size_t i, std::size_t j) const
{
	using std::acos;
	using std::fabs;

	T const result = fabs(phi(i) - phi(j));
	T const pi = acos(T(-1.0));

	return (result > pi) ? (T(2.0) * pi - result) : result;
}

template <typename T>
T phase_space_point<T>::cos_angle_neg(std::size_t i, std::size_t j) const
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
	T const cosh_rap = cosh(rapidity_shift_);
	T const sinh_rap = sinh(rapidity_shift_);
	T const pt = sqrt(px * px + py * py);
	T const qt = sqrt(qx * qx + qy * qy);

	T const pz = pt * (sinh_rap - cosh_rap * x1) / sqrt(T(1.0) - x1 * x1);
	T const qz = qt * (sinh_rap - cosh_rap * x2) / sqrt(T(1.0) - x2 * x2);

	T const num = px*qx + py*qy + pz*qz;
	T const den = sqrt((px*px + py*py + pz*pz) * (qx*qx + qy*qy + qz*qz));

	return num / den;
}

template <typename T>
T phase_space_point<T>::cos_angle_pos(std::size_t i, std::size_t j) const
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
	T const cosh_rap = cosh(rapidity_shift_);
	T const sinh_rap = sinh(rapidity_shift_);
	T const pt = sqrt(px * px + py * py);
	T const qt = sqrt(qx * qx + qy * qy);

	T const pz = pt * (sinh_rap + cosh_rap * x1) / sqrt(T(1.0) - x1 * x1);
	T const qz = qt * (sinh_rap + cosh_rap * x2) / sqrt(T(1.0) - x2 * x2);

	T const num = px*qx + py*qy + pz*qz;
	T const den = sqrt((px*px + py*py + pz*pz) * (qx*qx + qy*qy + qz*qz));

	return num / den;
}

template <typename T>
T phase_space_point<T>::dist2(std::size_t i, std::size_t j) const
{
	T const rap = rap_diff(i, j);
	T const phi = abs_phi_diff(i, j);

	return rap * rap + phi * phi;
}

template <typename T>
T phase_space_point<T>::m2(std::size_t i) const
{
	T const p0 = p_[4*i+0];
	T const p1 = p_[4*i+1];
	T const p2 = p_[4*i+2];
	T const p3 = p_[4*i+3];

	return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

template <typename T>
T phase_space_point<T>::m2(std::size_t i, std::size_t j) const
{
	T const p0 = p_[4*i+0] + p_[4*j+0];
	T const p1 = p_[4*i+1] + p_[4*j+1];
	T const p2 = p_[4*i+2] + p_[4*j+2];
	T const p3 = p_[4*i+3] + p_[4*j+3];

	return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

template <typename T>
T phase_space_point<T>::mt(std::size_t i, std::size_t j) const
{
	using std::sqrt;

	T const pt = sqrt(pt2(i)) + sqrt(pt2(j));
	T const px = p_[4*i+1] + p_[4*j+1];
	T const py = p_[4*i+2] + p_[4*j+2];

	return sqrt(pt * pt - px * px - py * py);
}

template <typename T>
T phase_space_point<T>::mt(
	std::size_t i,
	std::size_t j,
	std::size_t k,
	std::size_t l
) const {
	using std::sqrt;

	T const pt = sqrt(pt2(i))+sqrt(pt2(j))+sqrt(pt2(k))+sqrt(pt2(l));
	T const px = p_[4*i+1] + p_[4*j+1] + p_[4*k+1] + p_[4*l+1];
	T const py = p_[4*i+2] + p_[4*j+2] + p_[4*k+2] + p_[4*l+2];

	return sqrt(pt * pt - px * px - py * py);
}

template <typename T>
T phase_space_point<T>::phi(std::size_t i) const
{
	using std::atan2;

	T const p1 = p_[4*i+1];
	T const p2 = p_[4*i+2];

	return atan2(p2, p1);
}

template <typename T>
T phase_space_point<T>::pt2(std::size_t i) const
{
	T const p1 = p_[4*i+1];
	T const p2 = p_[4*i+2];

	return p1 * p1 + p2 * p2;
}

template <typename T>
T phase_space_point<T>::pt2(std::size_t i, std::size_t j) const
{
	T const p1 = p_[4*i+1] + p_[4*j+1];
	T const p2 = p_[4*i+2] + p_[4*j+2];

	return p1 * p1 + p2 * p2;
}

template <typename T>
T phase_space_point<T>::rap_diff(std::size_t i, std::size_t j) const
{
	using std::atanh;

	T const x = p_[4*i+3] / p_[4*i+0];
	T const y = p_[4*j+3] / p_[4*j+0];

	return atanh((x - y) / (T(1.0) - x * y));
}

template <typename T>
T phase_space_point<T>::rap_neg(std::size_t i) const
{
	using std::atanh;

	T const p0 = p_[4*i+0];
	T const p3 = p_[4*i+3];

	return rapidity_shift_ - atanh(p3 / p0);
}

template <typename T>
T phase_space_point<T>::rap_pos(std::size_t i) const
{
	using std::atanh;

	T const p0 = p_[4*i+0];
	T const p3 = p_[4*i+3];

	return rapidity_shift_ + atanh(p3 / p0);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class phase_space_point<double>;

}
