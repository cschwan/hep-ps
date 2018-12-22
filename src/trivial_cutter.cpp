#include "hep/ps/trivial_cutter.hpp"

#include <algorithm>

namespace hep
{

template <typename T>
bool trivial_cutter<T>::cut(psp<T> const&) const
{
    return false;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class trivial_cutter<double>;

}
