#include "hep/ps/phase_space_generator.hpp"

namespace hep
{

template <typename T>
phase_space_generator<T>::~phase_space_generator() = default;

template <typename T>
T phase_space_generator<T>::operator()(
	std::size_t channel,
	std::vector<T> const& random_numbers,
	std::vector<T>& momenta,
	std::vector<std::size_t> const&,
	std::vector<T>& densities,
	hep::multi_channel_map action
) {
	if (action == hep::multi_channel_map::calculate_densities)
	{
		return this->densities(densities);
	}

	generate(random_numbers, momenta, channel);

	// value gets ignored
	return T(1.0);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class phase_space_generator<double>;

}
