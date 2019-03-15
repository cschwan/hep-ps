#include "hep/ps/cs_subtraction.hpp"
#include "hep/ps/d_subtraction.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/phase_space_point.hpp"
#include "hep/ps/pdg_functions.hpp"

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

template <typename T>
void d_subtraction<T>::insertion_terms(
    int_dipole const& term,
    nonstd::span<scales<T> const> scales,
    std::vector<T> const& phase_space,
    T x,
    T eta,
    std::vector<ab_term<T>>& results
) const {
    subtraction_.insertion_terms(term, scales, phase_space, x, eta, results);

    auto const ex = pdg_id_to_particle_type(term.vertex().external());
    auto const in = pdg_id_to_particle_type(term.vertex().internal());

    if ((term.type() == insertion_term_type::born) && (ex == particle_type::fermion) &&
        (in == particle_type::fermion))
    {
        using std::acos;

        T const pi = acos(T(-1.0));

        for (std::size_t i = 0; i != scales.size(); ++i)
        {
            results.at(i).b -= T(0.5) / pi * (T(1.0) / T(3.0) * pi * pi - T(1.0));
        }
    }
}

template <typename T>
void d_subtraction<T>::insertion_terms2(
    int_dipole const& term,
    nonstd::span<scales<T> const> scales,
    std::vector<T> const& phase_space,
    std::vector<T>& results
) const {
    using std::acos;
    using std::log;

    results.clear();

    // TODO: the other cases are NYI
    assert( pdg_id_to_particle_type(term.vertex().unresolved()) == particle_type::boson );

    // TODO: QCD?
    assert( correction_type_of(term.vertex()) == correction_type::ew );

    T const sij = hep::phase_space_point<T>{phase_space}.m2(term.emitter(), term.spectator());
    T const pi = acos(T(-1.0));

    T dipole_dependent_factor = T();

    switch (term.type())
    {
    case insertion_term_type::initial_initial:
        dipole_dependent_factor -= pi * pi / T(3.0);
        dipole_dependent_factor += T(2.0);
        break;

    case insertion_term_type::initial_final:
        dipole_dependent_factor += pi * pi / T(6.0);
        dipole_dependent_factor -= T(1.0);
        break;

    case insertion_term_type::final_initial:
        dipole_dependent_factor -= pi * pi / T(2.0);
        dipole_dependent_factor += T(3.0) / T(2.0);
        break;

    case insertion_term_type::final_final:
        dipole_dependent_factor -= pi * pi / T(3.0);
        dipole_dependent_factor += T(3.0) / T(2.0);
        break;

    default:
        // implementation error
        assert( false );
    }

    for (auto const& scale : scales)
    {
        T const mu2 = scale.regularization() * scale.regularization();
        T const logmubsij = log(mu2 / sij);

        T result = T(2.0);
        result += T(0.5) * logmubsij * logmubsij;
        result += T(1.5) * logmubsij;
        result += dipole_dependent_factor;
        result *= T(-0.5) / pi;

        results.push_back(result);
    }
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class d_subtraction<double>;

}
