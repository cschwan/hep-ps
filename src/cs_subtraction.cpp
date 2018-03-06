#include "hep/ps/cs_subtraction.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/parton_type.hpp"
#include "hep/ps/phase_space_point.hpp"

#include <gsl/gsl_sf_dilog.h>

#include <cassert>
#include <cmath>
#include <cstddef>

namespace hep
{

template <typename T>
cs_subtraction<T>::cs_subtraction(
	T nc,
	T tf,
	T nf,
	factorization_scheme fscheme,
	regularization_scheme rscheme,
	insertion_term_mode mode
)
	: nc_{nc}
	, tf_{tf}
	, nf_{nf}
	, fscheme_{fscheme}
	, rscheme_{rscheme}
	, mode_{mode}
{
}

template <typename T>
dipole_invariants<T> cs_subtraction<T>::map_phase_space(
	std::vector<T> const& real_phase_space,
	std::vector<T>& born_phase_space,
	dipole const& dipole_info
) {
	std::size_t const i = dipole_info.emitter();
	std::size_t const j = dipole_info.unresolved();
	std::size_t const k = dipole_info.spectator();

	// 0. Assign the real momenta to the mapped phase space point, leaving out
	// the momentum of the unresolved particle.
	//

	born_phase_space.assign(
		real_phase_space.begin(),
		real_phase_space.begin() + 4 * j
	);
	born_phase_space.insert(
		born_phase_space.end(),
		real_phase_space.begin() + 4 * (j + 1),
		real_phase_space.end()
	);

	// 1. Calculate the invariants that are needed for both the phase space
	// mapping and the dipole functions.

	dipole_invariants<T> invariants;
	{
		phase_space_point<T> ps{real_phase_space};

		T const sij = ps.m2(i, j);
		T const sik = ps.m2(i, k);
		T const sjk = ps.m2(j, k);

		switch (dipole_info.type())
		{
		case dipole_type::final_final:
			invariants.one   = sij / (sij + sik + sjk);
			invariants.two   = sik / (      sik + sjk);
			invariants.alpha = invariants.one;

			break;

		case dipole_type::final_initial:
			invariants.alpha = sij / (sik + sjk);
			invariants.two   = sik / (sik + sjk);
			invariants.one   = T(1.0) - invariants.alpha;

			break;

		case dipole_type::initial_final:
			invariants.one   = T(1.0) - sjk / (sij + sik);
			invariants.two   =          sij / (sij + sik);
			invariants.alpha = invariants.two;

			break;

		case dipole_type::initial_initial:
			invariants.one   = T(1.0) - (sij + sjk) / sik;
			invariants.two   =           sij        / sik;
			invariants.alpha = invariants.two;

			break;
		}

		invariants.sij = sij;
	}

	// indices of recombined particle/spectator in the mapped phase space
	std::size_t ij_tilde = (j > i) ? i : i - 1;
	std::size_t k_tilde  = (j > k) ? k : k - 1;

	auto& pij0_tilde = born_phase_space.at(4 * ij_tilde + 0);
	auto& pij1_tilde = born_phase_space.at(4 * ij_tilde + 1);
	auto& pij2_tilde = born_phase_space.at(4 * ij_tilde + 2);
	auto& pij3_tilde = born_phase_space.at(4 * ij_tilde + 3);

	auto& pk0_tilde = born_phase_space.at(4 * k_tilde + 0);
	auto& pk1_tilde = born_phase_space.at(4 * k_tilde + 1);
	auto& pk2_tilde = born_phase_space.at(4 * k_tilde + 2);
	auto& pk3_tilde = born_phase_space.at(4 * k_tilde + 3);

	auto const& pj0 = real_phase_space.at(4 * j + 0);
	auto const& pj1 = real_phase_space.at(4 * j + 1);
	auto const& pj2 = real_phase_space.at(4 * j + 2);
	auto const& pj3 = real_phase_space.at(4 * j + 3);

	if (dipole_info.type() == dipole_type::final_final)
	{
		T const y = invariants.one;

		pk0_tilde /= T(1.0) - y;
		pk1_tilde /= T(1.0) - y;
		pk2_tilde /= T(1.0) - y;
		pk3_tilde /= T(1.0) - y;

		pij0_tilde += pj0 - y * pk0_tilde;
		pij1_tilde += pj1 - y * pk1_tilde;
		pij2_tilde += pj2 - y * pk2_tilde;
		pij3_tilde += pj3 - y * pk3_tilde;
	}
	else if (dipole_info.type() == dipole_type::final_initial)
	{
		T const x = invariants.one;

		pij0_tilde += pj0 - (T(1.0) - x) * pk0_tilde;
		pij1_tilde += pj1 - (T(1.0) - x) * pk1_tilde;
		pij2_tilde += pj2 - (T(1.0) - x) * pk2_tilde;
		pij3_tilde += pj3 - (T(1.0) - x) * pk3_tilde;

		pk0_tilde *= x;
		pk1_tilde *= x;
		pk2_tilde *= x;
		pk3_tilde *= x;
	}
	else if (dipole_info.type() == dipole_type::initial_final)
	{
		T const x = invariants.one;

		pk0_tilde += pj0 - (T(1.0) - x) * pij0_tilde;
		pk1_tilde += pj1 - (T(1.0) - x) * pij1_tilde;
		pk2_tilde += pj2 - (T(1.0) - x) * pij2_tilde;
		pk3_tilde += pj3 - (T(1.0) - x) * pij3_tilde;

		pij0_tilde *= x;
		pij1_tilde *= x;
		pij2_tilde *= x;
		pij3_tilde *= x;
	}
	else if (dipole_info.type() == dipole_type::initial_initial)
	{
		T const x = invariants.one;

		T const k0 = pij0_tilde + pk0_tilde - pj0;
		T const k1 = pij1_tilde + pk1_tilde - pj1;
		T const k2 = pij2_tilde + pk2_tilde - pj2;
		T const k3 = pij3_tilde + pk3_tilde - pj3;

		T const k0_tilde = x * pij0_tilde + pk0_tilde;
		T const k1_tilde = x * pij1_tilde + pk1_tilde;
		T const k2_tilde = x * pij2_tilde + pk2_tilde;
		T const k3_tilde = x * pij3_tilde + pk3_tilde;

		T const kk0 = k0 + k0_tilde;
		T const kk1 = k1 + k1_tilde;
		T const kk2 = k2 + k2_tilde;
		T const kk3 = k3 + k3_tilde;

		T const ksqr  = k0 * k0 - k1 * k1 - k2 * k2 - k3 * k3;
		T const kksqr = kk0 * kk0 - kk1 * kk1 - kk2 * kk2 - kk3 * kk3;

		for (std::size_t l = 0; l != born_phase_space.size() / 4; ++l)
		{
			if ((l == ij_tilde) || (l == k_tilde))
			{
				// skip the emitter momentum, we will set it at the end; also
				// skip the spectator momentum which remains unchanged
				continue;
			}

			T& p0 = born_phase_space.at(4 * l + 0);
			T& p1 = born_phase_space.at(4 * l + 1);
			T& p2 = born_phase_space.at(4 * l + 2);
			T& p3 = born_phase_space.at(4 * l + 3);

			T const pk  = p0 *  k0 - p1 *  k1 - p2 *  k2 - p3 *  k3;
			T const pkk = p0 * kk0 - p1 * kk1 - p2 * kk2 - p3 * kk3;

			p0 += T(2.0) * (pk / ksqr * k0_tilde - pkk / kksqr * kk0);
			p1 += T(2.0) * (pk / ksqr * k1_tilde - pkk / kksqr * kk1);
			p2 += T(2.0) * (pk / ksqr * k2_tilde - pkk / kksqr * kk2);
			p3 += T(2.0) * (pk / ksqr * k3_tilde - pkk / kksqr * kk3);
		}

		pij0_tilde *= x;
		pij1_tilde *= x;
		pij2_tilde *= x;
		pij3_tilde *= x;
	}
	else
	{
		// implementation error
		assert( false );
	}

	return invariants;
}

template <typename T>
T cs_subtraction<T>::fermion_function(
	dipole const& dipole_info,
	dipole_invariants<T> const& invariants
) {
	T const cf = tf_ * (nc_ * nc_ - T(1.0)) / nc_;
	T factor = T(8.0) * std::acos(T(-1.0)) * cf;

	T dipole;
	T propagator;

	switch (dipole_info.type())
	{
	case dipole_type::final_final:
	{
		T const y   = invariants.one;
		T const z   = invariants.two;
		T const sij = invariants.sij;

		dipole = T(2.0) / (T(1.0) - z * (T(1.0) - y)) - (T(1.0) + z);
		propagator = T(-1.0) / sij;
	}
		break;

	case dipole_type::final_initial:
	{
		T const x   = invariants.one;
		T const z   = invariants.two;
		T const sij = invariants.sij;

		dipole = T(2.0) / (T(1.0) - z + (T(1.0) - x)) - (T(1.0) + z);
		propagator = T(-1.0) / (sij * x);
	}
		break;

	case dipole_type::initial_final:
	{
		T const x   = invariants.one;
		T const u   = invariants.two;
		T const sij = invariants.sij;

		if (dipole_info.unresolved_type() == particle_type::fermion)
		{
			factor = T(8.0) * std::acos(T(-1.0)) * tf_;
			dipole = T(1.0) - T(2.0) * x * (T(1.0) - x);
		}
		else
		{
			dipole = T(2.0) / (T(1.0) - x + u) - (T(1.0) + x);
		}

		propagator = T(-1.0) / (sij * x);
	}
		break;

	case dipole_type::initial_initial:
	{
		T const x   = invariants.one;
		T const sij = invariants.sij;

		if (dipole_info.unresolved_type() == particle_type::fermion)
		{
			factor = T(8.0) * std::acos(T(-1.0)) * tf_;
			dipole = T(1.0) - T(2.0) * x * (T(1.0) - x);
		}
		else
		{
			dipole = T(2.0) / (T(1.0) - x) - (T(1.0) + x);
		}

		propagator = T(-1.0) / (sij * x);
	}
		break;

	default:
		// this signals an implementation error
		assert( false );
	}

	return factor * dipole * propagator;
}

template <typename T>
void cs_subtraction<T>::insertion_terms(
	insertion_term const& term,
	std::vector<scales<T>> const& scales,
	std::vector<T> const& phase_space,
	T x,
	T eta,
	std::vector<abc_terms<T>>& results
) const {
	using std::acos;
	using std::log;

	auto bo = parton_type::gluon_;

	switch (mode_)
	{
	case insertion_term_mode::ew:
	case insertion_term_mode::ew_no_photons:
	case insertion_term_mode::ew_only_photons:
		bo = parton_type::photon_;

		break;

	default:
		break;
	}

	constexpr auto qq = parton_type::quark;
	constexpr auto aq = parton_type::anti_quark;

	// TODO: DIS scheme is NYI
	assert( fscheme_ == factorization_scheme::msbar );

	T const cf = tf_ * (nc_ * nc_ - T(1.0)) / nc_;
	T const pi = acos(T(-1.0));

	results.clear();

	switch (term.type())
	{
	case insertion_term_type::born:
	{
		T const omx = T(1.0) - x;
		T const logomxbx = log(omx / x);
		T const dilogome = gsl_sf_dilog(T(1.0) - eta);
		T const logome = log(T(1.0) - eta);

		T const value2 = T(0.5) * cf / pi * (((T(1.0) + x * x) / omx) *
			logomxbx + omx);
		T const value3 = T(0.5) * cf / pi * (T(2.0) / omx) * logomxbx;
		T const value4 = T(0.5) * cf / pi * (T(2.0) / T(3.0) * pi * pi -
			T(5.0) + T(2.0) * dilogome + logome * logome);

		abc_terms<T> result;

		if (mode_ != insertion_term_mode::ew_no_photons)
		{
			T const value1 = T(0.5) * tf_ / pi * ((x * x + omx * omx) *
				logomxbx + T(2.0) * x * omx);

			result.a[bo][aq] = value1;
			result.a[bo][qq] = value1;
		}

		if (mode_ != insertion_term_mode::ew_only_photons)
		{
			result.a[aq][aq] = value2;
			result.a[qq][qq] = value2;
			result.b[aq][aq] = value3;
			result.b[qq][qq] = value3;
			result.c[aq][aq] = value4;
			result.c[qq][qq] = value4;
		}

		// TODO: qg and gg are NYI

		results.assign(scales.size(), result);
	}

		break;

	case insertion_term_type::final_initial:
	{
		T const ca = nc_;
		T const gamma = (term.emitter_type() == particle_type::fermion)
			? T(1.5) * cf
			: T(11.0) / T(6.0) * ca - T(2.0) / T(3.0) * tf_ * nf_;

		T const value1 = T(0.5) * gamma / pi / (T(1.0) - x);
		T const value2 = T(0.5) * gamma / pi * (T(1.0) + log(T(1.0) - eta));

		abc_terms<T> result;

		if (mode_ != insertion_term_mode::ew_only_photons)
		{
			result.a[aq][aq] = value1;
			result.a[qq][qq] = value1;
			result.b[aq][aq] = value1;
			result.b[qq][qq] = value1;
			result.c[aq][aq] = value2;
			result.c[qq][qq] = value2;
		}

		// TODO: qg and gg are NYI

		results.assign(scales.size(), result);
	}

		break;

	case insertion_term_type::initial_final:
	{
		phase_space_point<T> ps{phase_space};

		T const omx = T(1.0) - x;
		T const sai = ps.m2(term.emitter(), term.spectator());

		T const value1 = T(0.5) * tf_ / pi * (x * x + omx * omx);
		T const value2 = T(0.5) * cf / pi * (T(1.0) + x * x) / omx;
		T const value3 = T(0.5) * cf / pi * (T(0.5) * eta * (T(2.0) + eta) +
			T(2.0) * log(T(1.0) - eta));

		for (auto const& mu : scales)
		{
			T const mu2 = mu.factorization() * mu.factorization();
			T const logmu2bsai = log(mu2 / sai);

			abc_terms<T> result;

			if (mode_ != insertion_term_mode::ew_no_photons)
			{
				result.a[bo][aq] = value1 * logmu2bsai;
				result.a[bo][qq] = value1 * logmu2bsai;
			}

			if (mode_ != insertion_term_mode::ew_only_photons)
			{
				result.a[aq][aq] = value2 * logmu2bsai;
				result.a[qq][qq] = value2 * logmu2bsai;
				result.b[aq][aq] = value2 * logmu2bsai;
				result.b[qq][qq] = value2 * logmu2bsai;
				result.c[aq][aq] = value3 * logmu2bsai;
				result.c[qq][qq] = value3 * logmu2bsai;
			}

			// TODO: qg and gg are NYI

			results.push_back(result);
		}
	}

		break;

	case insertion_term_type::initial_initial:
	{
		phase_space_point<T> ps{phase_space};

		T const omx = T(1.0) - x;
		T const sai = ps.m2(term.emitter(), term.spectator());
		T const logomx = log(omx);
		T const logome = log(T(1.0) - eta);

		T const value1 = T(0.5) * tf_ / pi * (x * x + omx * omx);
		T const value2 = T(0.5) * cf / pi * (T(1.0) + x * x) / omx;
		T const value3 = T(0.5) * cf / pi * (T(0.5) * eta * (T(2.0) + eta) +
			T(2.0) * logome);
		T const value4 = T(0.5) * cf / pi * T(2.0) * logomx / omx;
		T const value5 = T(0.5) * cf / pi * (pi*pi / T(3.0) - logome * logome);

		for (auto const& mu : scales)
		{
			T const mu2 = mu.factorization() * mu.factorization();
			T const logmu2bsai = log(mu2 / sai);

			abc_terms<T> result;

			if (mode_ != insertion_term_mode::ew_no_photons)
			{
				result.a[bo][aq] = value1 * (logmu2bsai - logomx);
				result.a[bo][qq] = value1 * (logmu2bsai - logomx);
			}

			if (mode_ != insertion_term_mode::ew_only_photons)
			{
				result.a[aq][aq] = value2 * (logmu2bsai - logomx);
				result.a[qq][qq] = value2 * (logmu2bsai - logomx);
				result.b[aq][aq] = value2 * logmu2bsai - value4;
				result.b[qq][qq] = value2 * logmu2bsai - value4;
				result.c[aq][aq] = value3 * logmu2bsai + value5;
				result.c[qq][qq] = value3 * logmu2bsai + value5;
			}

			// TODO: qg and gg are NYI

			results.push_back(result);
		}

		break;
	}

	default:
		// implementation error
		assert( false );
	}
}

template <typename T>
void cs_subtraction<T>::insertion_terms2(
	insertion_term const& term,
	std::vector<scales<T>> const& scales,
	std::vector<T> const& phase_space,
	std::vector<T>& results
) const {
	using std::acos;

	if (term.type() == insertion_term_type::born)
	{
		results.assign(scales.size(), T());
		return;
	}

	T scheme_dep_factor;

	switch (rscheme_)
	{
	case regularization_scheme::dim_reg_blha:
		scheme_dep_factor = T(-1.0) / T(2.0);
		break;

	case regularization_scheme::dim_reg_coli:
		scheme_dep_factor = T(-2.0) / T(3.0);
		break;

	default:
		assert( false );
	}

	results.clear();

	switch (term.emitter_type())
	{
	case particle_type::fermion:
	{
		phase_space_point<T> ps{phase_space};

		T const pi = acos(T(-1.0));
		T const cf = tf_ * (nc_ * nc_ - T(1.0)) / nc_;
		T const sij = ps.m2(term.emitter(), term.spectator());

		for (auto const& mu : scales)
		{
			T const mu2 = mu.regularization() * mu.regularization();
			T const logmubsij = log(mu2 / sij);

			T result = T(5.0);
			result += scheme_dep_factor * pi * pi;
			result += T(1.5) * logmubsij;
			result += T(0.5) * logmubsij * logmubsij;
			result *= T(-0.5) * cf / pi;

			results.push_back(result);
		}
	}

		break;

	default:
		// NYI
		assert( false );
	}
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class cs_subtraction<double>;

}
