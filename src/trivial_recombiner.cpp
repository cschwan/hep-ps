#include "hep/ps/trivial_recombiner.hpp"

namespace hep
{

template <typename T>
std::size_t trivial_recombiner<T>::recombine(
	std::vector<T> const& phase_space,
	std::vector<T>& recombined_phase_space,
	std::vector<index_with_particle_class> const&,
	std::size_t
) const {
	recombined_phase_space = phase_space;

	return 0;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class trivial_recombiner<double>;

}
