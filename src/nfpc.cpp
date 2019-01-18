#include "hep/ps/ps_functions.hpp"
#include "hep/ps/nfpc.hpp"

#include "config.hpp"

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_dilog.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <utility>

namespace
{

enum invariant_type
{
    plain,
    tilde,
    bar,
    barbar
};

template <typename T>
T invariant(
    std::vector<T> const& phase_space,
    hep::nfpc_info<T> const& info,
    std::size_t one,
    std::size_t two,
    invariant_type type
) {
    T p0 = T();
    T p1 = T();
    T p2 = T();
    T p3 = T();

    switch (type)
    {
    case plain:
        p0 = phase_space.at(4 * one + 0) + phase_space.at(4 * two + 0);
        p1 = phase_space.at(4 * one + 1) + phase_space.at(4 * two + 1);
        p2 = phase_space.at(4 * one + 2) + phase_space.at(4 * two + 2);
        p3 = phase_space.at(4 * one + 3) + phase_space.at(4 * two + 3);

        break;

    case tilde:
        for (std::size_t const a : info.decays.at(one))
        {
            p0 += phase_space.at(4 * a + 0);
            p1 += phase_space.at(4 * a + 1);
            p2 += phase_space.at(4 * a + 2);
            p3 += phase_space.at(4 * a + 3);
        }

        p0 -= phase_space.at(4 * two + 0);
        p1 -= phase_space.at(4 * two + 1);
        p2 -= phase_space.at(4 * two + 2);
        p3 -= phase_space.at(4 * two + 3);

        break;

    case bar:
        for (std::size_t const a : info.decays.at(one))
        {
            p0 += phase_space.at(4 * a + 0);
            p1 += phase_space.at(4 * a + 1);
            p2 += phase_space.at(4 * a + 2);
            p3 += phase_space.at(4 * a + 3);
        }

        p0 += phase_space.at(4 * two + 0);
        p1 += phase_space.at(4 * two + 1);
        p2 += phase_space.at(4 * two + 2);
        p3 += phase_space.at(4 * two + 3);

        break;

    case barbar:
        for (std::size_t const a : info.decays.at(one))
        {
            p0 += phase_space.at(4 * a + 0);
            p1 += phase_space.at(4 * a + 1);
            p2 += phase_space.at(4 * a + 2);
            p3 += phase_space.at(4 * a + 3);
        }
        for (std::size_t const b : info.decays.at(two))
        {
            p0 += phase_space.at(4 * b + 0);
            p1 += phase_space.at(4 * b + 1);
            p2 += phase_space.at(4 * b + 2);
            p3 += phase_space.at(4 * b + 3);
        }

        break;

    default:
        assert( false );
    }

    T const result = p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3;

    return result;
}

template <typename T>
std::complex<T> off_shell_propagator(
    std::vector<T> const& phase_space,
    hep::nfpc_info<T> const& info,
    std::size_t i
) {
    T k0 = T();
    T k1 = T();
    T k2 = T();
    T k3 = T();

    for (std::size_t const a : info.resonance_decays.at(i))
    {
        k0 += phase_space.at(4 * a + 0);
        k1 += phase_space.at(4 * a + 1);
        k2 += phase_space.at(4 * a + 2);
        k3 += phase_space.at(4 * a + 3);
    }

    T const mass  = info.resonance_masses.at(i);
    T const width = info.resonance_widths.at(i);

    std::complex<T> const result(k0 * k0 - k1 * k1 - k2 * k2 - k3 * k3 - mass * mass, mass * width);

    return result;
}

template <typename T>
struct delta_invariants
{
    T mi;
    T mj;
    T mgamma;

    T sia;
    T sjb;
    T sib;
    T sja;
    T sij;
    T sab;

    std::complex<T> ki;
    std::complex<T> kj;
};

template <typename T>
std::complex<T> dilog(std::complex<T> const& argument)
{
    double const r     = static_cast <double> (std::abs(argument));
    double const theta = static_cast <double> (std::arg(argument));

    gsl_sf_result real;
    gsl_sf_result imag;

    auto const failed = gsl_sf_complex_dilog_e(r, theta, &real, &imag);

    assert( !failed );

    std::complex<T> result(T(real.val), T(imag.val));

    return result;
}

template <typename T>
std::complex<T> ac_dilog(std::complex<T> const& x1, std::complex<T> const& x2)
{
    using std::acos;
    using std::log;
    using std::signbit;

    auto const product = x1 * x2;

    bool const sign1 = signbit(x1.imag());
    bool const sign2 = signbit(x2.imag());
    bool const sign3 = signbit(product.imag());

    auto const result1 = dilog(T(1.0) - product);

    bool const pos =  sign1 &&  sign2 && !sign3;
    bool const neg = !sign1 && !sign2 &&  sign3;

    if (pos != neg)
    {
        auto const result2 = std::complex<T>(T(), (pos ? T(2.0) : T(-2.0)) * acos(T(-1.0))) *
            log(T(1.0) - product);

        return result1 + result2;
    }

    return result1;
}

// CHECKED!
template <typename T>
T delta_prime(delta_invariants<T> const& x)
{
    using std::log;

    std::complex<T> const ieps{ T(), T(1e-16) };

    T const mi2 = x.mi * x.mi;

    auto const arg1 = (mi2 - (x.sia + ieps)) / mi2;
    auto const arg2 = (mi2 - (x.sib + ieps)) / (-(x.sab + ieps));
    auto const arg3 = std::abs(-x.ki / (x.mgamma * x.mi));

    auto const dilog = ac_dilog(arg1, arg2);
    auto const result = (T(2.0) * (log(arg1) + log(arg2) - T(1.0)) * log(arg3) + T(2.0)) + dilog;

    return std::real(result);
}

// CHECKED!
template <typename T>
T delta_mm(delta_invariants<T> const& x)
{
    using std::log;

    auto const arg1 = (x.mgamma * x.mi) / -x.ki;
    auto const arg2 = (x.mgamma * x.mj) / -x.kj;
    auto const result = T(2.0) * log(arg1) + T(2.0) * log(arg2) + T(4.0);

    return std::real(result);
}

//// Implementation of Eq. (3.78) in Nucl. Phys. B 844 199-242
//template <typename T>
//std::complex<T> D00_new(delta_invariant<T> const& x)
//{
//  using std::sqrt;
//
//  T const m0 = x.mj;
//  T const m3 = x.mi;
//
//  T const Y01 = m0 * m0 - x.sjb;
//  T const Y02 = m0 * m0 - x.sja;
//  T const Y03 = m0 * m0 + m3 * m3 - x.sij;
//  T const Y12 = -x.sab;
//  T const Y13 = m3 * m3 - x.sib;
//  T const Y23 = m3 * m3 - x.sia;
//
//  T const a = Y13 * Y23 - m3 * m3 * Y12;
//  T const b = Y02 * Y13 + Y01 * Y23 - Y03 * Y12;
//  T const c = Y01 * Y02 - m0 * m0 * Y12;
//  T const d = Y12;
//
//  T const sqrt_detY = sqrt(b * b - T(4.0) * a * c);
//  T const x1 = (-b + sqrt_detY) * T(0.5) / a;
//  T const x2 = (-b - sqrt_detY) * T(0.5) / a;
//
//  T const root = sqrt(T(1.0) - T(4.0) * m0 * m3 /
//      (x.sij - (m0 - m3) * (m0 - m3)));
//  T const x03 = (root - T(1.0)) / (root + T(1.0));
//
//  T const r03_1 = m3 / m0 * x03;
//  T const r03_2 = m3 / m0 / x03;
//
//  std::complex<T> result;
//
//  result += ac_dilog(-x1, Y23 / Y02);
//  result -= ;
//  result -= ;
//  result += ;
//  result -= ;
//
//  result /= a * (x1 - x2);
//
//  return result;
//}

// Implementation of Eq. (C.3) in Nucl. Phys. B 519 39-84
template <typename T>
std::complex<T> D00(delta_invariants<T> const& x)
{
    using std::log;
    using std::sqrt;

    // TODO: this function is only valid for Mi == Mj and 1 -> 2 decays

    T const mi2 = x.mi * x.mi;
    T const mj2 = x.mj * x.mj;

    T const s13 = x.sib - x.sab - mi2;
    T const s14 = x.sij - x.sja - x.sib + x.sab;
    T const s23 = x.sab;
    T const s24 = x.sja - x.sab - mj2;

    std::complex<T> ieps{ T(), T(1e-16) };

    auto const kw = sqrt(hep::kaellen(mi2 * mi2, s13 * s24, s14 * s23) - ieps);

    auto const z = (mi2 * mi2 + s13 * s24 - s14 * s23 + kw) / (T(2.0) * s13 * s24);
    auto const x1 = s24 * z / mi2;
    auto const x2 = mi2 / (s13 * z);
    auto const bw = std::sqrt(T(1.0) - T(4.0) * mi2 / x.sij + ieps);
    auto const xw = (bw - T(1.0)) / (bw + T(1.0));

    std::complex<T> result;

    result -= ac_dilog(-(s13 + s23) / mi2 - ieps, -x1);
    result += ac_dilog(-(s13 + s23) / mi2 - ieps, -x2);
    result -= ac_dilog(-mi2 / (s23 + s24) + ieps, -x1);
    result += ac_dilog(-mi2 / (s23 + s24) + ieps, -x2);
    result += ac_dilog(xw, -x1);
    result -= ac_dilog(xw, -x2);
    result += ac_dilog(T(1.0) / xw, -x1);
    result -= ac_dilog(T(1.0) / xw, -x2);
    result += log(T(1.0) + s24 / s23) * log(-x1);
    result -= log(T(1.0) + s24 / s23) * log(-x2);

    result /= kw;

    return result;
}

template <typename T>
std::complex<T> D0_13(delta_invariants<T> const& x)
{
    using std::acos;
    using std::log;
    using std::sqrt;

    std::complex<T> ieps{ T(), T(1e-16) };

    auto const Ki  = x.ki;
    auto const Kj  = x.kj;
    auto const Mi  = x.mi;
    auto const Mi2 = x.mi * x.mi;
    auto const Mj  = x.mj;
    auto const Mj2 = x.mj * x.mj;
    auto const sab = x.sab + ieps;
    auto const sia = x.sia + ieps;
    auto const sib = x.sib + ieps;
    auto const sij = x.sij + ieps;
    auto const sjb = x.sjb + ieps;
    auto const ma  = (Mi2 - sia) / Mi;
    auto const mb  = (Mj2 - sjb) / Mj;

    assert( std::real(ma) > T() );
    assert( std::real(mb) > T() );

    // complex Kaellen function
    auto const sqrt_kaellen = sqrt((sij - Mj2 - Mi2) * (sij - Mj2 - Mi2) - T(4.0) * Mi2 * Mj2);

    auto const x1 = Ki / Kj;
    auto const x2 = (Mi2 - sib) / (Mj2 - sjb);
    auto const x3 = (Mi2 + Mj2 - sij + sqrt_kaellen) / (T(2.0) * Mi2);
    auto const x4 = (Mi2 + Mj2 - sij - sqrt_kaellen) / (T(2.0) * Mi2);

    std::complex<T> result;

    result -= T(2.0) * ac_dilog(x2, T(1.0) / x1);
    result += ac_dilog(T(1.0) / x1, T(1.0) / x3);
    result += ac_dilog(T(1.0) / x1, T(1.0) / x4);
    result -= ac_dilog(T(1.0) / x2, T(1.0) / x3);
    result -= ac_dilog(T(1.0) / x2, T(1.0) / x4);

    // arguments of the logarithms/dilogarithms
    auto const arg1 = (mb * Mi) / (Mi2 - sib);
    auto const arg3 = (Mi2 - sia) / Mi2;
    auto const arg4 = (Mi2 - sib) / (-sab);

    // actual logs
    auto const log1 = log(arg1);

    result += T(2.0) * log((ma * mb) / (-sab)) * log(x.mgamma * Mi / (-Ki));
    result -= log1 * log1;
    result -= ac_dilog(arg3, arg4);
    result -= acos(T(-1.0)) * acos(T(-1.0)) / T(3.0);

    return result / (sab * Ki);
}

template <typename T>
std::complex<T> D0_24(delta_invariants<T> const& x)
{
    using std::swap;

    delta_invariants<T> xprime(x);

    swap(xprime.mi, xprime.mj);
    swap(xprime.sia, xprime.sjb);
    swap(xprime.sib, xprime.sja);
    swap(xprime.ki, xprime.kj);

    return D0_13(xprime);
}

// Implementation of Eq. (C.1) in Nucl. Phys. B 519 39-84
template <typename T>
T delta_mm_prime(delta_invariants<T> const& x)
{
    using std::acos;
    using std::log;

    std::complex<T> ieps{ T() , T(1e-16) };

    auto const s  = x.sij;
    auto const mw = x.mi;
    auto const km = x.ki;
    auto const kp = x.kj;
    auto const lambda = x.mgamma;

    auto const bw = sqrt(T(1.0) - T(4.0) * mw * mw / s + ieps);
    auto const xw = (bw - T(1.0)) / (bw + T(1.0));

    T const pi = acos(T(-1.0));

    std::complex<T> result;

    result -= ac_dilog(km / kp,          xw);
    result += ac_dilog(km / kp, T(1.0) / xw);
    result +=    dilog(T(1.0) - xw * xw);
    result += pi * pi;
    result += log(-xw) * log(-xw);
    result += T(2.0) * log(kp / (lambda * mw)) * log(xw);
    result -= std::complex<T>{T(), T(2.0) * pi} * log(T(1.0) - xw * xw);

    result /= s * bw;

    T const factor = x.sij - x.mi * x.mi - x.mj * x.mj;

    return -factor * std::real(result);
}

// :s/Subscript(s,ab)/sab/ge
// :s/Subscript(K,i)/Ki/ge
// :s/Subscript(K,j)/Kj/ge
// :s/Power(Subscript(M,i),2)/Mi2/ge
// :s/Power(Subscript(M,j),2)/Mj2/ge
// :s/Subscript(OverBar(s),ja)/sja/ge
// :s/Subscript(OverBar(s),ib)/sib/ge
// :s/Subscript(OverBar(OverBar(s)),ij)/sij/ge
// :s/Power(Ki,2)/Ki*Ki/ge
// :s/Power(Kj,2)/Kj*Kj/ge
// :s/Power(Subscript(M,i),4)/Mi2*Mi2/ge
// :s/Power(sab,2)/sab*sab/ge
// :s/Power(\([^,]*\),\([0-9]*\)/pow(\1,T(\2.0)/ge
// :s/\(\<[0-9]\>\)\*/T(\1.0)\*/ge

template <typename T>
T delta_mf_prime_and_ff_prime(delta_invariants<T> const& x)
{
    // TODO: determinants are only valid for sia = sjb = 0

    using std::acos;
    using std::pow;

    auto const sab = x.sab;
    auto const sib = x.sib;
    auto const sja = x.sja;
    auto const sij = x.sij;
    auto const Mi2 = x.mi * x.mi;
    auto const Mj2 = x.mj * x.mj;
    auto const Ki  = x.ki;
    auto const Kj  = x.kj;

    auto const detY = T(2.0)*sab*(Kj*Kj*Mi2*(Mi2+sab-sib)+Ki*Ki*Mj2*(Mj2+sab-sja)+Ki*Kj*(Mj2*(-sab+
        sib)-sib*sja+Mi2*(-T(2.0)*Mj2-sab+sja)+sab*sij));

    auto const detY0 = -(Mi2*Mi2*pow(sab-sja,T(2.0)))-pow(Mj2*(sab-sib)+sib*sja-sab*sij,T(2.0))+
        T(2.0)*Mi2*(Mj2*(sab*sab-sib*sja-sab*(sib+sja-T(2.0)*sij))+(sab-sja)*(-(sib*sja)+sab*sij));

    auto const detY2 = sab*(T(2.0)*Ki*Mj2*(Mj2+sab-sja)+Kj*(Mj2*(-sab+sib)-sib*sja+Mi2*(-T(2.0)*Mj2-
        sab+sja)+sab*sij));

    auto const detY3 = sab*(T(2.0)*Kj*Mi2*(Mi2+sab-sib)+Ki*(Mj2*(-sab+sib)-sib*sja+Mi2*(-T(2.0)*Mj2-
        sab+sja)+sab*sij));

    auto const factor = sab * Ki * Kj;
    auto const result0 = detY0 * D00(x);
    auto const result2 = detY2 * D0_24(x);
    auto const result3 = detY3 * D0_13(x);

    auto const nominator = result0 + result2 + result3;
    auto const result    = nominator / detY;

    T const pi = acos(T(-1.0));

    return -std::real(factor * result + pi * pi / T(3.0));
}

}

namespace hep
{

template <typename T>
T nfpc(
    std::vector<T> const& phase_space,
    std::vector<T> const& resonance_invariants,
    nfpc_info<T> const& info
) {
    std::size_t const r = info.decays.size();

    T result = T();

    delta_invariants<T> x;
    x.mgamma = info.photon_mass;

    for (std::size_t i = 0; i != r; ++i)
    {
    for (auto const a : info.decays.at(i))
    {
        x.mi = info.masses.at(i);
        x.sia = invariant(phase_space, info, i, a, tilde);
        x.ki = std::complex<T>(resonance_invariants.at(i) - x.mi * x.mi, x.mi * info.widths.at(i));

        T const signa = info.signs.at(a);
        T const qa = info.charges.at(a);

        if (qa == T())
        {
            continue;
        }

        for (std::size_t j = i + 1; j != r; ++j)
        {
            x.mj = info.masses.at(j);
            x.sja = invariant(phase_space, info, j, a, bar);
            x.sij = invariant(phase_space, info, i, j, barbar);
            x.kj = std::complex<T>(resonance_invariants.at(j) - x.mj * x.mj, x.mj *
                info.widths.at(i));

            for (auto const b : info.decays.at(j))
            {
                T const signb = info.signs.at(b);
                T const qb = info.charges.at(b);

                if (qb == T())
                {
                    continue;
                }

                x.sjb = invariant(phase_space, info, j, b, tilde);
                x.sib = invariant(phase_space, info, i, b, bar);
                x.sab = invariant(phase_space, info, a, b, plain);

                // mm
                T delta1 = T();
                // mm'
                T delta2 = T();
                // mf + mf' + ff'
                T delta3 = T();

                delta1 += delta_mm(x);
                delta2 += delta_mm_prime(x);
                delta3 += delta_mf_prime_and_ff_prime(x);

                result += signa * signb * qa * qb * (delta1 + delta2 + delta3);
            }
        }

        for (auto const b : info.non_resonant_particles)
        {
            T const signb = info.signs.at(b);
            T const qb = info.charges.at(b);

            if (qb == T())
            {
                continue;
            }

            x.sib = invariant(phase_space, info, i, b, bar);
            x.sab = invariant(phase_space, info, a, b, plain);

            T const delta = delta_prime(x);

            result += signa * signb * qa * qb * delta;
        }
    }
    }

    return -result;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template double nfpc(
    std::vector<double> const&,
    std::vector<double> const&,
    nfpc_info<double> const&
);

}
