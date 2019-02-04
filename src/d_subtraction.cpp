#include "hep/ps/cs_subtraction.hpp"
#include "hep/ps/d_subtraction.hpp"

#include <cmath>

namespace hep
{

template <typename T>
d_subtraction<T>::d_subtraction(
    T nc,
    T tf,
    std::size_t nf,
    factorization_scheme fscheme,
    regularization_scheme rscheme
)
    : subtraction_{nc, tf, nf, fscheme, rscheme}
{
}

template <typename T>
bool d_subtraction<T>::same_mapping(dipole const& a, dipole const& b) const
{
    return subtraction_.same_mapping(a, b);
}

template <typename T>
dipole_invariants<T> d_subtraction<T>::map_phase_space(
    std::vector<T> const& real_phase_space,
    std::vector<T>& born_phase_space,
    dipole const& dipole_info
) {
    return subtraction_.map_phase_space(real_phase_space, born_phase_space, dipole_info);
}

template <typename T>
spin_correlation_matrix<T> d_subtraction<T>::boson_function(
    dipole const& dipole_info,
    dipole_invariants<T> const& invariants,
    std::vector<T> const& phase_space
) {
    return subtraction_.boson_function(dipole_info, invariants, phase_space);
}

template <typename T>
T d_subtraction<T>::fermion_function(
    dipole const& dipole_info,
    dipole_invariants<T> const& invariants
) {
    // TODO: not checked are dipoles with initial photons

    T result = subtraction_.fermion_function(dipole_info, invariants);

    if (dipole_info.type() == dipole_type::final_final)
    {
        // Stefan's FF dipoles have an additional 1/(1-y) factor
        result /= T(1.0) - invariants.one;
    }

    // FI, IF, and II dipoles are exactly the same

    return result;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class d_subtraction<double>;

}
