#include "hep/ps/p_type_jet_algorithm.hpp"
#include "hep/ps/phase_space_point.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>

namespace hep
{

template <typename T>
p_type_jet_algorithm<T>::p_type_jet_algorithm(T p, T radius)
	: p_{p}
	, radius2_{radius * radius}
{
}

template <typename T>
bool p_type_jet_algorithm<T>::find_jet(std::vector<T>& phase_space)
{
	using std::fmin;
	using std::pow;

	std::size_t const n = candidates_.size();

	if (n == 1)
	{
		return true;
	}

	dib_.clear();

	phase_space_point<T> ps{phase_space};

	for (std::size_t i = 0; i != n; ++i)
	{
		std::size_t const index = candidates_.at(i);
		dib_.push_back(pow(ps.pt2(index), p_));
	}

	auto const min_dib = std::min_element(dib_.begin(), dib_.end());
	std::size_t const min_ib = std::distance(dib_.begin(), min_dib);

	T min_dij = std::numeric_limits<T>::max();
	std::size_t min_i = 0;
	std::size_t min_j = 0;

	for (std::size_t i = 0; i != (n-1); ++i)
	{
		T const dib_ib = dib_.at(i);
		std::size_t const index_i = candidates_.at(i);

		for (std::size_t j = i + 1; j < n; ++j)
		{
			T const dib_jb = dib_.at(j);
			T const min = fmin(dib_ib, dib_jb);
			std::size_t const index_j = candidates_.at(j);
			T const factor = ps.dist2(index_i, index_j) / radius2_;
			T const dij = min * factor;

			if (dij < min_dij)
			{
				min_dij = dij;
				min_i = i;
				min_j = j;
			}
		}
	}

	if (*min_dib < min_dij)
	{
		candidates_.erase(std::next(candidates_.begin(), min_ib));

		return true;
	}

	auto const index_i = candidates_.at(min_i);
	auto const index_j = candidates_.at(min_j);

	phase_space.at(4 * index_i + 0) += phase_space.at(4 * index_j + 0);
	phase_space.at(4 * index_i + 1) += phase_space.at(4 * index_j + 1);
	phase_space.at(4 * index_i + 2) += phase_space.at(4 * index_j + 2);
	phase_space.at(4 * index_i + 3) += phase_space.at(4 * index_j + 3);

	phase_space.erase(
		std::next(phase_space.begin(), 4 * index_j + 0),
		std::next(phase_space.begin(), 4 * index_j + 4)
	);

	auto begin = std::next(candidates_.begin(), min_j);
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
	std::vector<final_state> const& final_states,
	std::size_t max_recombinations
) {
	std::size_t recombinations = 0;

	// we are only interested in jets, so we need to pick out the partons
	candidates_.clear();
	for (std::size_t i = 0; i != final_states.size(); ++i)
	{
		if (final_states.at(i) == final_state::quark_gluon)
		{
			candidates_.push_back(i + 2);
		}
	}

	recombined_phase_space = phase_space;

	dib_.reserve(candidates_.size());

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
