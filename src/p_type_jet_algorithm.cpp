#include "hep/ps/p_type_jet_algorithm.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>

namespace
{

template <typename T>
T rap(std::vector<T> const& p, std::size_t i)
{
	T const p0 = p[4 * i + 0];
	T const p3 = p[4 * i + 3];

	return T(0.5) * std::log((p0 + p3) / (p0 - p3));
}

template <typename T>
T phi(std::vector<T> const& p, std::size_t i)
{
	T const p1 = p[4 * i + 1];
	T const p2 = p[4 * i + 2];

	return std::atan2(p2, p1);
}

template <typename T>
T phi_diff(T phi1, T phi2)
{
	T const abs_delta_phi = std::fabs(phi1 - phi2);

	if (abs_delta_phi > std::acos(T(-1.0)))
	{
		return T(2.0) * std::acos(T(-1.0)) - abs_delta_phi;
	}

	return abs_delta_phi;
}

template <typename T>
T dist2(std::vector<T> const& p, std::size_t i, std::size_t j)
{
	T const rap_diff = rap(p, i) - rap(p, j);
	T const angle_diff = phi_diff(phi(p, i), phi(p, j));

	return rap_diff * rap_diff + angle_diff * angle_diff;
}

template <typename T>
T pt2(std::vector<T> const& p, std::size_t i)
{
	T const p1 = p.at(4 * i + 1);
	T const p2 = p.at(4 * i + 2);

	return p1 * p1 + p2 * p2;
}

constexpr std::ptrdiff_t reconstruct_distance(
	std::size_t n,
	std::size_t i,
	std::size_t j
) {
	return ((2 * n - i - 3) * i) / 2 + j - 1;
}

}

namespace hep
{

template <typename T>
p_type_jet_algorithm<T>::p_type_jet_algorithm(T p, T radius)
	: p_(p)
	, radius_(radius)
{
}

template <typename T>
bool p_type_jet_algorithm<T>::find_jet(std::vector<T>& phase_space)
{
	std::size_t const n = candidates_.size();

	if (n == 1)
	{
		return true;
	}

	dib_.clear();
	dij_.clear();

	for (std::size_t i = 0; i != n; ++i)
	{
		std::size_t const index = candidates_.at(i);
		dib_.push_back(std::pow(pt2(phase_space, index), p_));
	}

	for (std::size_t i = 0; i != (n-1); ++i)
	{
		T const dib_ib = dib_.at(i);
		std::size_t const index_i = candidates_.at(i);

		for (std::size_t j = i + 1; j < n; ++j)
		{
			T const dib_jb = dib_.at(j);
			T const min = std::fmin(dib_ib, dib_jb);
			std::size_t const index_j = candidates_.at(j);
			T const factor = dist2(phase_space, index_i, index_j) /
				(radius_ * radius_);

			dij_.push_back(min * factor);
		}
	}

	auto const min_dib = std::min_element(dib_.begin(), dib_.end());
	auto const min_dij = std::min_element(dij_.begin(), dij_.end());

	if (*min_dib < *min_dij)
	{
		std::size_t const index = std::distance(dib_.begin(), min_dib);
		candidates_.erase(std::next(candidates_.begin(), index));

		return true;
	}

	auto const distance = std::distance(dij_.begin(), min_dij);

	std::size_t i = 0;
	std::size_t j = 1;

	// TODO: find a way to remove the loop
	while ((distance - reconstruct_distance(n, i, j)) >=
		static_cast <std::ptrdiff_t> (n - i - 1))
	{
		++i;
		++j;
	}

	j += distance - reconstruct_distance(n, i, j);

	auto const index_i = candidates_.at(i);
	auto const index_j = candidates_.at(j);

	phase_space.at(4 * index_i + 0) += phase_space.at(4 * index_j + 0);
	phase_space.at(4 * index_i + 1) += phase_space.at(4 * index_j + 1);
	phase_space.at(4 * index_i + 2) += phase_space.at(4 * index_j + 2);
	phase_space.at(4 * index_i + 3) += phase_space.at(4 * index_j + 3);

	phase_space.erase(
		std::next(phase_space.begin(), 4 * index_j + 0),
		std::next(phase_space.begin(), 4 * index_j + 4)
	);

	auto begin = std::next(candidates_.begin(), j);
	auto next = std::next(begin);
	auto end = candidates_.end();

	std::transform(next, end, next, [](std::size_t v) { return v - 1; });
	std::rotate(begin, next, end);
	candidates_.erase(std::prev(end));

	return false;
}

template <typename T>
std::size_t p_type_jet_algorithm<T>::recombine(
	std::vector<T> const& phase_space,
	std::vector<T>& recombined_phase_space,
	std::vector<std::size_t> const& recombination_candidates,
	std::size_t max_recombinations
) {
	std::size_t recombinations = 0;
	candidates_ = recombination_candidates;

	recombined_phase_space = phase_space;

	std::size_t const n = candidates_.size();
	dib_.reserve(n);
	dij_.reserve((n * (n - 1)) >> 1);

	for (;;)
	{
		if (find_jet(recombined_phase_space))
		{
			if (candidates_.size() == 1)
			{
				break;
			}
		}
		else
		{
			++recombinations;

			if (recombinations > max_recombinations)
			{
				break;
			}
		}
	}

	return recombinations;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class p_type_jet_algorithm<double>;

}
