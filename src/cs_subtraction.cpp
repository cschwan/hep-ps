#include "hep/ps/cs_subtraction.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

namespace
{

template <typename T>
T invariant(std::vector<T> const& p, std::size_t i, std::size_t j)
{
	T const p0 = p.at(4 * i + 0) + p.at(4 * j + 0);
	T const p1 = p.at(4 * i + 1) + p.at(4 * j + 1);
	T const p2 = p.at(4 * i + 2) + p.at(4 * j + 2);
	T const p3 = p.at(4 * i + 3) + p.at(4 * j + 3);

	return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;
}

}

namespace hep
{

template <typename T>
cs_subtraction<T>::cs_subtraction(T n, T tf)
	: n_(n)
	, tf_(tf)
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
		T const sij = invariant(real_phase_space, i, j);
		T const sik = invariant(real_phase_space, i, k);
		T const sjk = invariant(real_phase_space, j, k);

		switch (dipole_info.type())
		{
		case dipole_type::final_final:
			invariants.one     = sij / (sij + sik + sjk);
			invariants.two     = sik / (      sik + sjk);
			invariants.adipole = invariants.one;

			break;

		case dipole_type::final_initial:
			invariants.adipole = sij / (sik + sjk);
			invariants.two     = sik / (sik + sjk);
			invariants.one     = T(1.0) - invariants.adipole;

			break;

		case dipole_type::initial_final:
			invariants.one     = T(1.0) - sjk / (sij + sik);
			invariants.two     =          sij / (sij + sik);
			invariants.adipole = invariants.two;

			break;

		case dipole_type::initial_initial:
			invariants.one     = T(1.0) - (sij + sjk) / sik;
			invariants.two     =           sij        / sik;
			invariants.adipole = invariants.two;

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
	T const cf = tf_ * (n_ * n_ - T(1.0)) / n_;
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

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class cs_subtraction<double>;

}