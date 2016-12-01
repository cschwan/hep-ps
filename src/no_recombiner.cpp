#include "hep/ps/no_recombiner.hpp"

namespace hep
{

template <typename T>
std::size_t no_recombiner<T>::recombine(
	std::vector<T> const& phase_space,
	std::vector<T>& recombined_phase_space,
	std::vector<std::size_t> const&,
	std::size_t
) const {
	recombined_phase_space = phase_space;

	return 0;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class no_recombiner<double>;

}
