#include "hep/ps/no_cutter.hpp"

namespace hep
{

template <typename T>
cut_result no_cutter<T>::cut(
	std::vector<T> const&,
	T,
	bool
) {
	// cut nothing
	return cut_result(false, false);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class no_cutter<double>;

}
