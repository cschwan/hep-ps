#include "hep/ps/cs_subtraction.hpp"
#include "hep/ps/d_subtraction.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/phase_space_point.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <gsl/gsl_sf_dilog.h>

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
    using std::acos;
    using std::log;

    // TODO: the other cases are NYI
    assert( pdg_id_to_particle_type(term.vertex().unresolved()) == particle_type::boson );

    // TODO: QCD NYI and probably not needed
    assert( correction_type_of(term.vertex()) == correction_type::ew );

    T const pi = acos(T(-1.0));

    results.clear();

    switch (term.type())
    {
    case insertion_term_type::born:
    {
        T const omx = T(1.0) - x;
        T const opxs = T(1.0) + x * x;
        T const logomx = log(omx);
        T const logome = log(T(1.0) - eta);
        T const s = phase_space_point<T>{phase_space}.m2(0, 1);

        ab_term<T> result = {};

        for (auto const& mu : scales)
        {
            T const muf = mu.factorization();
            T const logmuf2bs = log(muf * muf / s);

            T const value1 = opxs / omx * (logmuf2bs - T(2.0) * logomx + T(1.0));
            T const value2 = T(2.0) * eta - logome * (eta * (T(2.0) + eta) - T(1.0)
                + T(2.0) * logome) + T(0.5) * (eta * (T(2.0) + eta) + T(4.0) * logome) * logmuf2bs;

            result.a = T(-0.5) / pi *  value1;
            result.b = T(-0.5) / pi * (value1 + (eta - T(1.0)) * value2);

            results.push_back(result);
        }
    }

        break;

    case insertion_term_type::final_initial:
    {
        T const omx = T(1.0) - x;
        T const logomx = log(omx);
        T const dilogemo = gsl_sf_dilog(eta - T(1.0));
        T const logome = log(T(1.0) - eta);

        ab_term<T> result = {};

        T const value1 = T(1.0) / omx * (T(2.0) * log(T(2.0) - x) - T(2.0) * logomx
            - T(3.0) / T(2.0));
        T const value2 = -pi * pi / T(6.0) - T(3.0) / T(2.0) * logome - logome * logome
            - T(2.0) * dilogemo;

        result.a = T(-0.5) / pi * value1;
        result.b = result.a + (eta - T(1.0)) *T(-0.5) / pi *  value2;

        results.assign(scales.size(), result);
    }

        break;

    case insertion_term_type::initial_final:
    {
        T const omx = T(1.0) - x;
        T const pff = (T(1.0) + x * x) / omx;
        T const logomx = log(omx);
        T const logtmx = log(T(2.0) - x);
        T const s = phase_space_point<T>{phase_space}.m2(0, 1);
        T const sia = phase_space_point<T>{phase_space}.m2(term.emitter(), term.spectator());

        ab_term<T> result = {};

        T const value1 = pff * (log(sia / (s * x)) - T(1.0)) - T(2.0) / omx * logtmx
            + (T(1.0) + x) * logomx + omx;
        T const value2 = T();

        result.a = T(-0.5) / pi * value1;
        result.b = result.a + (eta - T(1.0)) * T(-0.5) / pi * value2;

        results.assign(scales.size(), result);
    }

        break;

    case insertion_term_type::initial_initial:
    {
        T const omx = T(1.0) - x;
        T const pff = (T(1.0) + x * x) / omx;
        T const s = phase_space_point<T>{phase_space}.m2(0, 1);
        T const sia = phase_space_point<T>{phase_space}.m2(term.emitter(), term.spectator());

        T const dilogome = gsl_sf_dilog(T(1.0) - eta);
        T const logome = log(T(1.0) - eta);

        ab_term<T> result = {};

        T const value1 = pff * (log(T(1.0) / x) - T(1.0)) - omx;
        T const value2 = pi * pi / T(3.0) - T(2.0) * dilogome - eta * (T(0.25) * eta
            + T(0.5) * (T(2.0) + eta) * log(eta)) + T(2.0) * logome - log(sia / s)
            * (T(2.0) * logome + eta + T(0.5) * eta * eta);

        result.a = T(-0.5) / pi * value1;
        result.b = result.a + (eta - T(1.0)) * T(-0.5) / pi * value2;

        results.assign(scales.size(), result);
    }
        break;

    default:
        assert( false );
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
