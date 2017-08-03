#include "hep/ps/phase_space_generator.hpp"

namespace hep
{

template <typename T>
phase_space_generator<T>::~phase_space_generator() = default;

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class phase_space_generator<double>;

}
