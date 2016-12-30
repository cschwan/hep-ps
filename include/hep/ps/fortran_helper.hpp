#ifndef HEP_PS_FORTRAN_HELPER_HPP
#define HEP_PS_FORTRAN_HELPER_HPP

#include <vector>

namespace hep
{

/// Assumes `momenta` contains a phase space point in FORTRAN ordering, i.e.
/// first all the energies, then all x-components, and so on. Reorders `momenta`
/// to C++ ordering, i.e. where the momenta are in `0 1 2 3`, `0 1 2 3`, ...
/// ordering.
template <typename T>
void fortran_ordering_to_cpp(std::vector<T>& momenta);

/// Has the opposite effect of \ref fortran_ordering_to_cpp.
template <typename T>
void cpp_ordering_to_fortran(std::vector<T>& momenta);

}

#endif
