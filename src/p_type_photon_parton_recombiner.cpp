#include "hep/ps/p_type_photon_parton_recombiner.hpp"
#include "hep/ps/phase_space_point.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace hep
{

template <typename T>
p_type_photon_parton_recombiner<T>::p_type_photon_parton_recombiner(
	T p,
	T photon_radius,
	T parton_radius
)
	: p_{p}
	, photon_radius2_{photon_radius * photon_radius}
	, jet_algorithm_{p, parton_radius}
{
}

template <typename T>
void p_type_photon_parton_recombiner<T>::recombine(
	std::vector<T> const& phase_space,
	std::vector<final_state> const& final_states,
	std::vector<T>& recombined_phase_space,
	std::vector<recombined_state>& recombined_states
) {
	using std::fmin;
	using std::pow;

	std::size_t const photons = std::count(
		final_states.begin(),
		final_states.end(),
		final_state::photon
	);

	// TODO: implement the general case
	assert( photons < 2 );

	fermions_.clear();
	photons_.clear();

	std::vector<final_state> final_state_copy(final_states);

	for (std::size_t i = 0; i != final_states.size(); ++i)
	{
		switch (final_states.at(i))
		{
		case final_state::charged_lepton:
			fermions_.push_back(i + 2);
			break;

		case final_state::quark_gluon:
			fermions_.push_back(i + 2);
			break;

		case final_state::photon:
			photons_.push_back(i + 2);
			break;

		case final_state::neutrino:
			break;
		}
	}

	recombined_phase_space = phase_space;

	dib_fermions_.reserve(fermions_.size());
	dib_photons_.reserve(photons_.size());

	while ((photons > 0) && (fermions_.size() > 0))
	{
		dib_photons_.clear();
		dib_fermions_.clear();

		phase_space_point<T> ps{recombined_phase_space};

		for (auto const fermion : fermions_)
		{
			dib_fermions_.push_back(pow(ps.pt2(fermion), p_));
		}
		for (auto const photon : photons_)
		{
			dib_photons_.push_back(pow(ps.pt2(photon), p_));
		}

		auto const min_dib_fermion = std::min_element(
			dib_fermions_.begin(),
			dib_fermions_.end()
		);
		std::size_t const min_ib_fermion = std::distance(
			dib_fermions_.begin(),
			min_dib_fermion
		);
		auto const min_dib_photon = std::min_element(
			dib_photons_.begin(),
			dib_photons_.end()
		);

		T min_dij = std::numeric_limits<T>::max();
		std::size_t min_i = 0;
		std::size_t min_j = 0;

		for (std::size_t i = 0; i != fermions_.size(); ++i)
		{
			T const dib_ib = dib_fermions_.at(i);
			std::size_t const index_i = fermions_.at(i);

			for (std::size_t j = 0; j != photons_.size(); ++j)
			{
				T const dib_jb = dib_photons_.at(j);
				T const min = fmin(dib_ib, dib_jb);
				std::size_t const index_j = photons_.at(j);
				T const factor = ps.dist2(index_i, index_j) / photon_radius2_;
				T const dij = min * factor;

				if (dij < min_dij)
				{
					min_dij = dij;
					min_i = i;
					min_j = j;
				}
			}
		}

		if ((*min_dib_photon < min_dij) || (*min_dib_fermion < min_dij))
		{
			if (*min_dib_photon < *min_dib_fermion)
			{
				break;
			}
			else
			{
				fermions_.erase(std::next(fermions_.begin(), min_ib_fermion));
			}
		}
		else
		{
			auto const index_i = fermions_.at(min_i);
			auto const index_j = photons_.at(min_j);

			recombined_phase_space.at(4 * index_i + 0) +=
				recombined_phase_space.at(4 * index_j + 0);
			recombined_phase_space.at(4 * index_i + 1) +=
				recombined_phase_space.at(4 * index_j + 1);
			recombined_phase_space.at(4 * index_i + 2) +=
				recombined_phase_space.at(4 * index_j + 2);
			recombined_phase_space.at(4 * index_i + 3) +=
				recombined_phase_space.at(4 * index_j + 3);

			recombined_phase_space.erase(
				std::next(recombined_phase_space.begin(), 4 * index_j + 0),
				std::next(recombined_phase_space.begin(), 4 * index_j + 4)
			);

			final_state_copy.erase(std::next(final_state_copy.begin(),
				index_j - 2));

			// early exit because we are only interested in a single photon
			break;

			// TODO: implement the general case

//			auto begin = std::next(candidates_.begin(), min_j);
//			auto next = std::next(begin);
//			auto end = candidates_.end();
//
//			std::transform(next, end, next, [](std::size_t p) { return p- 1; });
//			std::rotate(begin, next, end);
//			candidates_.erase(std::prev(end));
		}
	}

	jet_algorithm_.recombine(
		std::vector<T>(recombined_phase_space),
		final_state_copy,
		recombined_phase_space,
		recombined_states
	);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class p_type_photon_parton_recombiner<double>;

}
