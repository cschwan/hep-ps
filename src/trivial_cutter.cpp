#include "hep/ps/trivial_cutter.hpp"

#include <algorithm>

namespace hep
{

template <typename T>
cut_result trivial_cutter<T>::cut(std::vector<T> const&, T, std::vector<recombined_state> const&)
{
    // cut nothing
    return cut_result(false, false);
}

template <typename T>
bool trivial_cutter<T>::cut(psp<T> const&) const
{
    return false;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class trivial_cutter<double>;

}
