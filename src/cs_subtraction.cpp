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
	std::size_t nf,
	factorization_scheme fscheme,
	regularization_scheme rscheme
)
	: nc_{nc}
	, tf_{tf}
	, nf_{nf}
	, fscheme_{fscheme}
	, rscheme_{rscheme}
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
spin_correlation_matrix<T> cs_subtraction<T>::boson_function(
	dipole const& dipole_info,
	dipole_invariants<T> const& invariants,
	std::vector<T> const& phase_space
) {
	using std::acos;

	spin_correlation_matrix<T> result;

	std::size_t const i = dipole_info.emitter();
	std::size_t const j = dipole_info.unresolved();
	std::size_t const k = dipole_info.spectator();

	T const factor = T(8.0) * acos(T(-1.0));

	// FIXME: color is probably missing in the QCD case

	// TODO: gluon splitting into gluons NYI
	assert( (dipole_info.emitter_type() == particle_type::fermion) &&
		(dipole_info.unresolved_type() == particle_type::fermion) );

	switch (dipole_info.type())
	{
	case dipole_type::final_final:
	case dipole_type::final_initial:
	{
		T const zi = invariants.two;
		T const zj = T(1.0) - zi;

		result.a = factor * T(-1.0) / invariants.sij;
		result.b = factor * T(-1.0) / invariants.sij * T(4.0) * zi * zj;

		if (dipole_info.type() == dipole_type::final_initial)
		{
			result.a /= invariants.one;
			result.b /= invariants.one;
		}

		result.p = {
			phase_space.at(4 * i + 0) * zi - phase_space.at(4 * j + 0) * zj,
			phase_space.at(4 * i + 1) * zi - phase_space.at(4 * j + 1) * zj,
			phase_space.at(4 * i + 2) * zi - phase_space.at(4 * j + 2) * zj,
			phase_space.at(4 * i + 3) * zi - phase_space.at(4 * j + 3) * zj,
		};

	}

		break;

	case dipole_type::initial_final:
	case dipole_type::initial_initial:
	{
		T const x = invariants.one;

		result.a = factor * T(-1.0) / invariants.sij;
		result.b = factor * T(-1.0) / invariants.sij * T(4.0) * (x - T(1.0)) /
			(x * x);

		T const uj = invariants.two;
		T const uk = (dipole_info.type() == dipole_type::initial_initial)
			? T(1.0) : (T(1.0) - uj);

		result.p = {
			phase_space.at(4 * j + 0) * uk - phase_space.at(4 * k + 0) * uj,
			phase_space.at(4 * j + 1) * uk - phase_space.at(4 * k + 1) * uj,
			phase_space.at(4 * j + 2) * uk - phase_space.at(4 * k + 2) * uj,
			phase_space.at(4 * j + 3) * uk - phase_space.at(4 * k + 3) * uj,
		};

	}

		break;

	default:
		assert( false );
	}


	return result;
}

template <typename T>
T cs_subtraction<T>::fermion_function(
	dipole const& dipole_info,
	dipole_invariants<T> const& invariants
) {
	using std::acos;

	T const cf = tf_ * (nc_ * nc_ - T(1.0)) / nc_;
	T factor = T(8.0) * acos(T(-1.0));

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
			// TODO: check that the `if` was missing previously
			if (dipole_info.corr_type() == correction_type::qcd)
			{
				factor *= tf_ / cf;
			}

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
			if (dipole_info.corr_type() == correction_type::qcd)
			{
				factor *= tf_ / cf;
			}

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
	std::vector<ab_terms<T>>& results
) const {
	using std::acos;
	using std::log;

	parton_type const bo = (term.corr_type() == correction_type::ew)
		? parton_type::photon_
		: parton_type::gluon_;

	constexpr auto qq = parton_type::quark;
	constexpr auto aq = parton_type::anti_quark;

	// TODO: DIS scheme is NYI
	assert( fscheme_ == factorization_scheme::msbar );

	T const pi = acos(T(-1.0));
	T const color = (term.corr_type() == correction_type::ew)
		? nc_
		: (nc_ / (nc_ * nc_ - T(1.0))); // tf/cf

	T const color2 = (term.corr_type() == correction_type::ew)
		? T(1.0)
		: tf_ * (T(1.0) - T(1.0) / nc_ / nc_); // cf/ca

	results.clear();

	switch (term.type())
	{
	case insertion_term_type::born:
	{
		T const omx = T(1.0) - x;
		T const logomxbx = log(omx / x);
		T const dilogome = gsl_sf_dilog(T(1.0) - eta);
		T const logome = log(T(1.0) - eta);

		T const value1 = T(0.5) * color / pi * ((x * x + omx * omx) *
			logomxbx + T(2.0) * x * omx);
		T const value2 = T(0.5) / pi * (((T(1.0) + x * x) / omx) * logomxbx +
			omx);
		T const value3 = T(0.5) / pi * (T(2.0) / omx) * logomxbx;
		T const value4 = T(0.5) / pi * (T(2.0) / T(3.0) * pi * pi - T(5.0) +
			T(2.0) * dilogome + logome * logome);
		T const value5 = T(0.5) * color2 / pi * ((T(1.0) + omx * omx) / x *
			logomxbx + x);

		ab_terms<T> result;

		result.a[bo][aq] = value1;
		result.a[bo][qq] = value1;
		result.a[aq][aq] = value2;
		result.a[qq][qq] = value2;
		result.b[aq][aq] = (eta - T(1.0)) * value3 + value4;
		result.b[qq][qq] = (eta - T(1.0)) * value3 + value4;
		result.a[aq][bo] = value5;
		result.a[qq][bo] = value5;

		// TODO: gg is NYI

		results.assign(scales.size(), result);
	}

		break;

	case insertion_term_type::final_initial:
	{
		T const ca = nc_;
		T const gamma = (term.emitter_type() == particle_type::fermion)
			? T(1.5)
			// FIXME: the following branch has probably the wrong color factors
			: T(11.0) / T(6.0) * ca - T(2.0) / T(3.0) * tf_ * T(nf_);

		T const value1 = T(0.5) * gamma / pi / (T(1.0) - x);
		T const value2 = T(0.5) * gamma / pi * (T(1.0) + log(T(1.0) - eta));

		ab_terms<T> result;

		result.a[aq][aq] = value1;
		result.a[qq][qq] = value1;
		result.b[aq][aq] = (eta - T(1.0)) * value1 + value2;
		result.b[qq][qq] = (eta - T(1.0)) * value1 + value2;

		// TODO: gg is NYI

		results.assign(scales.size(), result);
	}

		break;

	case insertion_term_type::initial_final:
	{
		phase_space_point<T> ps{phase_space};

		T const omx = T(1.0) - x;
		T const sai = ps.m2(term.emitter(), term.spectator());

		T const value1 = T(0.5) * color / pi * (x * x + omx * omx);
		T const value2 = T(0.5) / pi * (T(1.0) + x * x) / omx;
		T const value3 = T(0.5) / pi * (T(0.5) * eta * (T(2.0) + eta) +
			T(2.0) * log(T(1.0) - eta));
		T const value4 = T(0.5) * color2 / pi * (T(1.0) + omx * omx) / x;

		for (auto const& mu : scales)
		{
			T const mu2 = mu.factorization() * mu.factorization();
			T const logmu2bsai = log(mu2 / sai);

			T const result1 = logmu2bsai * value1;
			T const result2 = logmu2bsai * value2;
			T const result3 = logmu2bsai * ((eta - T(1.0)) * value2 + value3);
			T const result4 = logmu2bsai * value4;

			ab_terms<T> result;

			result.a[bo][aq] = result1;
			result.a[bo][qq] = result1;
			result.a[aq][aq] = result2;
			result.a[qq][qq] = result2;
			result.b[aq][aq] = result3;
			result.b[qq][qq] = result3;
			result.a[aq][bo] = result4;
			result.a[qq][bo] = result4;

			// TODO: gg is NYI

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

		T const value1 = T(0.5) * color / pi * (x * x + omx * omx);
		T const value2 = T(0.5) / pi * (T(1.0) + x * x) / omx;
		T const value3 = T(0.5) / pi * (T(0.5) * eta * (T(2.0) + eta) +
			T(2.0) * logome);
		T const value4 = T(0.5) / pi * T(2.0) * logomx / omx;
		T const value5 = T(0.5) / pi * (pi*pi / T(3.0) - logome * logome);
		T const value6 = T(0.5) * color2 / pi * (T(1.0) + omx * omx) / x;

		for (auto const& mu : scales)
		{
			T const mu2 = mu.factorization() * mu.factorization();
			T const logmu2bsai = log(mu2 / sai);

			T const result1 = value1 * (logmu2bsai - logomx);
			T const result2 = value2 * (logmu2bsai - logomx);
			T const result3 = (eta - T(1.0)) * (value2 * logmu2bsai - value4) +
				value3 * logmu2bsai + value5;
			T const result4 = value6 * (logmu2bsai - logomx);

			ab_terms<T> result;

			result.a[bo][aq] = result1;
			result.a[bo][qq] = result1;
			result.a[aq][aq] = result2;
			result.a[qq][qq] = result2;
			result.b[aq][aq] = result3;
			result.b[qq][qq] = result3;
			result.a[aq][bo] = result4;
			result.a[qq][bo] = result4;

			// TODO: gg is NYI

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
		// should not be called with this function
		assert( false );
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

	phase_space_point<T> ps{phase_space};

	T const pi = acos(T(-1.0));
	T const sij = ps.m2(term.emitter(), term.spectator());

	switch (term.emitter_type())
	{
	case particle_type::boson:
	{
		if (term.corr_type() == correction_type::ew)
		{
			std::array<T, 6> const ncq2 = {
				// up
				nc_ * (T(0.0) * T(1.0) / T(9.0) + T(1.0) * T(4.0) / T(9.0)),
				// up, down
				nc_ * (T(1.0) * T(1.0) / T(9.0) + T(1.0) * T(4.0) / T(9.0)),
				// up, down, strange
				nc_ * (T(2.0) * T(1.0) / T(9.0) + T(1.0) * T(4.0) / T(9.0)),
				// up, down, strange, charm
				nc_ * (T(2.0) * T(1.0) / T(9.0) + T(2.0) * T(4.0) / T(9.0)),
				// up, down, strange, charm, bottom
				nc_ * (T(3.0) * T(1.0) / T(9.0) + T(2.0) * T(4.0) / T(9.0)),
				// up, down, strange, charm, bottom, top
				nc_ * (T(3.0) * T(1.0) / T(9.0) + T(3.0) * T(4.0) / T(9.0))
			};

			T const gamma = T(-2.0) / T(3.0) * ncq2.at(nf_);

			for (auto const& mu : scales)
			{
				T const mu2 = mu.regularization() * mu.regularization();
				T const logmubsij = log(mu2 / sij);

				// for photons there is no 1/eps^2 pole -> BHLA/COLI are equal
				T result = T(8.0) / T(3.0) * gamma;
				result += gamma * logmubsij;
				// TODO: shouldn't this be -1/2?
				result *= T(0.5) / pi;

				results.push_back(result);
			}
		}
		else if (term.corr_type() == correction_type::qcd)
		{
			T const ca = nc_;
			T const trnfbca = tf_ * nf_ / ca;

			for (auto const& mu : scales)
			{
				T const mu2 = mu.regularization() * mu.regularization();
				T const logmubsij = log(mu2 / sij);

				T result = T(223.0) / T(18.0) - T(16.0) / T(9.0) * trnfbca;
				result += scheme_dep_factor * pi * pi;
				result -= T(2.0) / T(3.0) * trnfbca * logmubsij;
				result += T(0.5) * logmubsij * logmubsij;
				result *= T(-0.5) / pi;

				results.push_back(result);
			}
		}
		else
		{
			assert( false );
		}
	}

		break;

	case particle_type::fermion:
	{
		for (auto const& mu : scales)
		{
			T const mu2 = mu.regularization() * mu.regularization();
			T const logmubsij = log(mu2 / sij);

			T result = T(5.0);
			result += scheme_dep_factor * pi * pi;
			result += T(1.5) * logmubsij;
			result += T(0.5) * logmubsij * logmubsij;
			result *= T(-0.5) / pi;

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
