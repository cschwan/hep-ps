#include "hep/ps/trivial_cutter.hpp"

namespace hep
{

template <typename T>
cut_result trivial_cutter<T>::cut(std::vector<T> const&, T, bool)
{
	// cut nothing
	return cut_result(false, false);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class trivial_cutter<double>;

}
