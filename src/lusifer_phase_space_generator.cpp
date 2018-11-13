#include "hep/ps/boost.hpp"
#include "hep/ps/fortran_helper.hpp"
#include "hep/ps/kaellen.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "hep/ps/lusifer_ps_channels.hpp"

#include "hadron_hadron_psg_adapter.hpp"
#include "lusifer_interfaces.hpp"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>
#include <utility>

namespace
{

constexpr bool included(int binary1, int binary2)
{
    return (binary1 & binary2) == binary2;
}

struct invariant
{
    invariant(std::size_t in, std::size_t idhep)
        : in(in)
        , idhep(idhep)
    {
    }

    std::size_t in;
    std::size_t idhep;
    std::bitset<std::numeric_limits<std::size_t>::digits> lmin;
    std::bitset<std::numeric_limits<std::size_t>::digits> lmax;
    std::size_t index;
};

struct process
{
    process(
        std::size_t in1,
        std::size_t in2,
        std::size_t out1,
        std::size_t out2,
        std::size_t in,
        std::size_t virt,
        std::size_t idhep
    )
        : in1(in1)
        , in2(in2)
        , out1(out1)
        , out2(out2)
        , in(in)
        , virt(virt)
        , idhep(idhep)
    {
    }

    std::size_t in1;
    std::size_t in2;
    std::size_t out1;
    std::size_t out2;
    std::size_t in;
    std::size_t virt;
    std::size_t idhep;
    std::size_t index;
};

struct decay
{
    decay(
        std::size_t in,
        std::size_t out1,
        std::size_t out2
    )
        : in(in)
        , out1(out1)
        , out2(out2)
    {
    }

    std::size_t in;
    std::size_t out1;
    std::size_t out2;
    std::size_t index;
};

bool operator==(decay const& a, decay const& b)
{
    if (a.in != b.in)
    {
        return false;
    }

    if ((a.out1 != b.out1) && (a.out1 != b.out2))
    {
        return false;
    }

    return true;
}

struct channel
{
    std::vector<invariant> invariants;
    std::vector<process> processes;
    std::vector<decay> decays;
};

struct invariant_info
{
    invariant_info(std::size_t index, std::size_t channel)
        : index(index)
        , channel(channel)
    {
    }

    std::size_t index;
    std::size_t channel;
};

struct process_info
{
    process_info(std::size_t index, std::size_t channel)
        : index(index)
        , channel(channel)
    {
    }

    std::size_t index;
    std::size_t channel;
};

struct decay_info
{
    decay_info(std::size_t index, std::size_t channel)
        : index(index)
        , channel(channel)
    {
    }

    std::size_t index;
    std::size_t channel;
};

template <typename T>
struct particle_info
{
    particle_info(T power, T mass, T width)
        : power(power)
        , mass(mass)
        , width(width)
    {
    }

    T power;
    T mass;
    T width;
};

template <typename T>
T jacobian(T power, T mass, T width, T x, T xmin, T xmax)
{
    using std::atan;
    using std::copysign;
    using std::fabs;
    using std::log;
    using std::pow;

    T m2 = copysign(mass * mass, mass);
    mass = fabs(mass);

    if (width > T())
    {
        T const mg = mass * width;
        T const min = (xmin - m2) / mg;
        T const max = (xmax - m2) / mg;

        // formula for the difference of two arctans
        T max_min = atan((xmax - xmin) / mg / (T(1.0) + max * min));
        if (min * max < T(-1.0))
        {
            max_min += (max > T()) ? acos(T(-1.0)) : -acos(T(-1.0));
        }

        T const xprime = x - m2;

        return mg / (max_min * (xprime * xprime + mg * mg));
    }

    if (power == T())
    {
        return T(1.0) / (xmax - xmin);
    }

    if (mass == T())
    {
        m2 -= T(1e-6);
    }

    if (power == T(1.0))
    {
        // TODO: WARNING this branch is untested!
        return T(1.0) / (log((xmax - m2) / (xmin - m2)) * (x - m2));
    }

    T const omp = T(1.0) - power;

    return omp / ((pow(fabs(xmax - m2), omp) - pow(fabs(xmin - m2), omp)) *
        pow(fabs(x - m2), power));
}

template <typename T>
T map(T power, T mass, T width, T x, T xmin, T xmax)
{
    using std::atan;
    using std::copysign;
    using std::exp;
    using std::fabs;
    using std::log;
    using std::pow;
    using std::tan;

    T m2 = copysign(mass * mass, mass);
    mass = fabs(mass);
    T result;

    if (width > T())
    {
        T const mg = mass * width;
        T const min = (xmin - m2) / mg;
        T const max = (xmax - m2) / mg;

        // formula for the difference of two arctans
        T max_min = atan((xmax - xmin) / mg / (T(1.0) + max * min));

        if (min * max < T(-1.0))
        {
            max_min += (max > T()) ? acos(T(-1.0)) : -acos(T(-1.0));
        }

        result = m2 + mg * tan(x * max_min + atan(min));
    }
    else if (power == T())
    {
        result = x * xmax + (T(1.0) - x) * xmin;

        // capture this case and return to avoid infinite recursion
        if ((result < xmin) || (result > xmax))
        {
            // returning the upper bound produces the least problems
            return xmax;
        }
    }
    else
    {
        if (mass == T())
        {
            // this reduces the number of extremely small invariants which
            // introduce numerical problems in the matrix elements
            m2 -= T(1e-6);
        }

        if (power == T(1.0))
        {
            // TODO: WARNING this branch is untested!
            result = m2 + exp(x * log(xmax - m2) + (T(1.0) - x) *
                log(xmin - m2));
        }
        else
        {
            T const omp = T(1.0) - power;

            result = m2 + pow(x * pow(fabs(xmax - m2), omp) + (T(1.0) - x) *
                pow(fabs(xmin - m2), omp), T(1.0) / omp);
        }
    }

    // if the result is outside of the expected interval, calculate it linearly;
    // this will only happen if `xmax` is ridiculously small or if `xmin` is
    // close to `xmax`, in which case it is a valid approximation
    if ((result < xmin) || (result > xmax))
    {
        return map(T{}, T{}, T{}, x, xmin, xmax);
    }

    return result;
}

template <typename T>
void rotate(std::array<T, 4>& p, T phi, T cos_theta)
{
    using std::cos;
    using std::fabs;
    using std::sin;
    using std::sqrt;

    T const sin_phi = sin(phi);
    T const cos_phi = cos(phi);
    T const sin_theta = sqrt(fabs((T(1.0) - cos_theta) * (T(1.0) + cos_theta)));

    T const px = p[1];
    T const py = p[2];
    T const pz = p[3];

    p[1] =  px * cos_theta * cos_phi + py * sin_phi + pz * sin_theta * cos_phi;
    p[2] = -px * cos_theta * sin_phi + py * cos_phi - pz * sin_theta * sin_phi;
    p[3] = -px * sin_theta                          + pz * cos_theta;
}

template <typename T>
void decay_momenta(
    T m,
    std::array<T, 4> const& q,
    T phi,
    T cos_theta,
    std::array<T, 4>& p1,
    std::array<T, 4>& p2
) {
    rotate(p1, phi, cos_theta);
    hep::boost(m, p1, q, true);

    p2[0] = q[0] - p1[0];
    p2[1] = q[1] - p1[1];
    p2[2] = q[2] - p1[2];
    p2[3] = q[3] - p1[3];
}

template <typename T>
struct tinv
{
    T tmin;
    T tmax;
    T lambdas;
    T lambdat;
};

template <typename T>
tinv<T> calc_tinv(T s, T s1, T s2, T t1, T t2)
{
    using std::fabs;
    using std::sqrt;

    T const threshold = T(1e-5);

    std::size_t const non_zero_invariants =
        ((s1 == T{}) ? 0 : 1) | ((s2 == T{}) ? 0 : 2) |
        ((t1 == T{}) ? 0 : 4) | ((t2 == T{}) ? 0 : 8);

    T const lambdas = hep::sqrt_kaellen(s, s1, s2);
    T const lambdat = hep::sqrt_kaellen(s, t1, t2);

    T tmin = T{};
    T tmax = T{};

    T const x = s - s1 - s2;
    T const y = s - t1 - t2;
    T const e1 = s1 * s2 / (x * x);
    T const e2 = t1 * t2 / (y * y);

    switch (non_zero_invariants)
    {
    case 1: // s2 = t1 = t2 = 0
    case 2: // s1 = t1 = t2 = 0
        tmin = -lambdas;
        // tmax is zero
        break;

    case 3: // t1 = t2 = 0
        tmin = T(-0.5) * x * (T(1.0) + sqrt(fabs(T(1.0) - T(4.0) * e1)));
        tmax = T(-0.5) * x * (T(1.0) - sqrt(fabs(T(1.0) - T(4.0) * e1)));
        // TODO: separate code path for small e1
        break;

    case 4: // s1 = s2 = t2 = 0
    case 8: // s1 = s2 = t1 = 0
        tmin = -lambdat;
        // tmax is zero
        break;

    case 5: // s2 = t2 = 0
        tmin = -s + s1 + t1 - s1 * t1 / s;
        // tmax is zero
        break;

    case 6: // s1 = t2 = 0
        tmin = -s + s2 + t1;
        tmax = s2 * t1 / s;
        break;

    case 9: // s2 = t1 = 0
        tmin = -s + s1 + t2;
        tmax = s1 * t2 / s;
        break;

    case 10: // s1 = t1 = 0
        tmin = -s + s2 + t2 - s2 * t2 / s;
        // tmax is zero
        break;

    case 12: // s1 = s2 = 0
        tmin = T(-0.5) * y * (T(1.0) + sqrt(fabs(T(1.0) - T(4.0) * e2)));
        tmax = T(-0.5) * y * (T(1.0) - sqrt(fabs(T(1.0) - T(4.0) * e2)));
        // TODO: separate code path for small e2
        break;

    default:
        T tmp = (s + s1 - s2) * (s + t1 - t2);
        tmin = s1 + t1 - T(0.5) * (tmp + lambdas * lambdat) / s;
        tmax = s1 + t1 - T(0.5) * (tmp - lambdas * lambdat) / s;

        if ((e1 < threshold) && (e2 < threshold))
        {
            // TODO: add order six terms?
            tmax = ((t2 * s1 + t1 * s2) - (x * y * (e1 + e2 + (e1 - e2) * (e1 - e2)))) / s;
        }
    }

    return { tmin, tmax, lambdas, lambdat };
}

template <typename T>
class lusifer_psg
{
public:
    using numeric_type = T;

    lusifer_psg(
        std::vector<std::string> const& processes,
        hep::lusifer_constants<T> const& constants,
        std::size_t extra_random_numbers = 0
    );

    std::size_t channels() const;

    T densities(std::vector<T>& densities);

    std::size_t dimensions() const;

    void generate(
        std::vector<T> const& random_numbers,
        std::vector<T>& momenta,
        T cmf_energy,
        std::size_t channel
    );

    std::size_t map_dimensions() const;

private:
    bool invariants_equal(int ch1, int ns1, int ch2, int ns2, int nex);
    bool processes_equal(int ch1, int ns1, int ch2, int ns2, int nex);

    std::vector<channel> channels_;
    std::vector<invariant_info> invariants;
    std::vector<process_info> processes_;
    std::vector<decay_info> decays;
    std::vector<particle_info<T>> particle_infos;
    std::vector<T> mcut;
    std::vector<T> invariant_jacobians;
    std::vector<T> process_jacobians;
    std::vector<T> decay_jacobians;
    std::vector<T> s;
    std::vector<std::array<T, 4>> p;
    std::size_t particles;
    std::size_t extra_random_numbers;
    T cmf_energy_;
};

template <typename T>
bool lusifer_psg<T>::invariants_equal(int ch1, int ns1, int ch2, int ns2, int nex)
{
    auto const& a = channels_.at(ch1).invariants.at(ns1);
    auto const& b = channels_.at(ch2).invariants.at(ns2);

    if (a.in != b.in)
    {
        return false;
    }

    auto const i = particle_infos.at(a.idhep);
    auto const j = particle_infos.at(b.idhep);

    if (i.power != j.power)
    {
        return false;
    }

    if (i.mass != j.mass)
    {
        return false;
    }

    // TODO: the following statement is probably unnecessary
    if (i.width != j.width)
    {
        return false;
    }

    for (int i = 0; i < (nex - 5); ++i)
    {
        if (a.lmin.test(i))
        {
            bool result = false;

            for (int j = 0; j < (nex - 5); ++j)
            {
                if (b.lmin.test(j) && (channels_.at(ch1).invariants.at(i).in ==
                    channels_.at(ch2).invariants.at(j).in))
                {
                    result = true;
                }
            }

            if (!result)
            {
                return false;
            }
        }

        if (a.lmax.test(i))
        {
            bool result = false;

            for (int j = 0; j < (nex - 5); ++j)
            {
                if (b.lmax.test(j) && (channels_.at(ch1).invariants.at(i).in ==
                    channels_.at(ch2).invariants.at(j).in))
                {
                    result = true;
                }
            }

            if (!result)
            {
                return false;
            }
        }
    }

    for (int j = 0; j < (nex - 5); ++j)
    {
        if (b.lmin.test(j))
        {
            bool result = false;

            for (int i = 0; i < (nex - 5); ++i)
            {
                if (a.lmin.test(i) && (channels_.at(ch1).invariants.at(i).in ==
                    channels_.at(ch2).invariants.at(j).in))
                {
                    result = true;
                }
            }

            if (!result)
            {
                return false;
            }
        }

        if (b.lmax.test(j))
        {
            bool result = false;

            for (int i = 0; i < (nex - 5); ++i)
            {
                if (a.lmax.test(i) && (channels_.at(ch1).invariants.at(i).in ==
                    channels_.at(ch2).invariants.at(j).in))
                {
                    result = true;
                }
            }

            if (!result)
            {
                return false;
            }
        }
    }

    // ignore member `index`, since it is calculated using this function

    return true;
}

template <typename T>
bool lusifer_psg<T>::processes_equal(int ch1, int ns1, int ch2, int ns2, int nex)
{
    auto const& a = channels_.at(ch1).processes.at(ns1);
    auto const& b = channels_.at(ch2).processes.at(ns2);

    if (a.in != b.in)
    {
        return false;
    }

    if ((a.virt != b.virt) && (a.virt != ((1 << nex) - 3 - b.virt)))
    {
        return false;
    }

    auto const i = particle_infos.at(a.idhep);
    auto const j = particle_infos.at(b.idhep);

    if (i.power != j.power)
    {
        return false;
    }

    if (i.mass != j.mass)
    {
        return false;
    }

    // TODO: the following statement is probably unnecessary
    if (i.width != j.width)
    {
        return false;
    }

    return true;
}

template <typename T>
lusifer_psg<T>::lusifer_psg(
    std::vector<std::string> const& processes,
    hep::lusifer_constants<T> const& constants,
    std::size_t extra_random_numbers
) {
    auto const result = hep::lusifer_ps_channels(processes, constants);
    auto const nex = processes.at(0).size() / 3;

    s.resize(1 << nex);
    p.resize(1 << nex);

    this->particles = nex;
    this->extra_random_numbers = extra_random_numbers;
    this->channels_.resize(result.size());

    for (std::size_t i = 0; i != result.size(); ++i)
    {
        auto& channel = this->channels_.at(i);

        channel.invariants.reserve(result.at(i).invariants().size());
        channel.processes.reserve(result.at(i).tchannels().size());
        channel.decays.reserve(result.at(i).decays().size());

        for (std::size_t j = 0; j != channel.invariants.capacity(); ++j)
        {
            channel.invariants.emplace_back(
                result.at(i).invariants().at(j).mom_id(),
                result.at(i).invariants().at(j).pdg_id()
            );
        }

        for (std::size_t a = 0; a != channel.invariants.size(); ++a)
        {
            for (std::size_t b = 0; b != a; ++b)
            {
                std::size_t const binary1 = channel.invariants.at(a).in + 1;
                std::size_t const binary2 = channel.invariants.at(b).in + 1;

                bool bit = false;

                if (included(binary1, binary2))
                {
                    bit = true;

                    for (std::size_t c = 0; c != a; ++c)
                    {
                        std::size_t const binary3 = channel.invariants.at(c).in + 1;

                        if (included(binary3, binary2) && (b != c))
                        {
                            bit = false;
                        }
                    }
                }

                assert( bit == lusifer_cinv.lmin[0][i][a][b] );

                channel.invariants.at(a).lmin.set(b, bit);
            }
        }

        std::size_t const allbinary = (1 << nex) - 1;

        for (std::size_t a = 0; a != channel.invariants.size(); ++a)
        {
            for (std::size_t b = 0; b != a; ++b)
            {
                std::size_t const binary1 = allbinary - 3 - channel.invariants.at(a).in - 1;
                std::size_t const binary2 = channel.invariants.at(b).in + 1;

                bool bit = false;

                if (included(binary1, binary2))
                {
                    bit = true;

                    for (std::size_t c = 0; c != a; ++c)
                    {
                        std::size_t binary3 = channel.invariants.at(c).in + 1;

                        if (included(binary3, binary2) && (b != c))
                        {
                            bit = false;
                        }
                    }
                }

                assert( bit == lusifer_cinv.lmax[0][i][a][b] );

                channel.invariants.at(a).lmax.set(b, bit);
            }
        }

        for (std::size_t j = 0; j != channel.processes.capacity(); ++j)
        {
            channel.processes.emplace_back(
                result.at(i).tchannels().at(j).in1_mom_id(),
                result.at(i).tchannels().at(j).in2_mom_id(),
                result.at(i).tchannels().at(j).out1_mom_id(),
                result.at(i).tchannels().at(j).out2_mom_id(),
                result.at(i).tchannels().at(j).in_mom_id(),
                result.at(i).tchannels().at(j).virt_id(),
                result.at(i).tchannels().at(j).pdg_id()
            );
        }

        for (std::size_t j = 0; j != channel.decays.capacity(); ++j)
        {
            channel.decays.emplace_back(
                result.at(i).decays().at(j).in_mom_id(),
                result.at(i).decays().at(j).out1_mom_id(),
                result.at(i).decays().at(j).out2_mom_id()
            );
        }
    }

    particle_infos.reserve(maxv + 1);

    for (std::size_t i = 0; i != particle_infos.capacity(); ++i)
    {
        // TODO: make this parameter available from outside
        T power = T(0.9);
        T mass  = lusifer_general.mass[i];
        T width = lusifer_general.width[i];

        if ((width > T()) || (i == 0) || ((i >= 30) && (i < 33)))
        {
            power = T();
        }

        particle_infos.emplace_back(power, mass, width);
    }

    std::size_t counter = 0;

    for (std::size_t a = 0; a != channels_.size(); ++a)
    {
        for (std::size_t b = 0; b != channels_.at(a).invariants.size(); ++b)
        {
            std::size_t index = 0;

            for (std::size_t c = 0; c != a; ++c)
            {
                for (std::size_t d = 0; d != channels_.at(c).invariants.size(); ++d)
                {
                    if (invariants_equal(a, b, c, d, nex))
                    {
                        index = channels_.at(c).invariants.at(d).index + 1;

                        // break out of two loops
                        c = a - 1;
                        break;
                    }
                }
            }

            if (index == 0)
            {
                ++counter;
                index = counter;

                assert( a == (lusifer_cdensity.chinv[0][invariants.size()] - 1) );
                assert( b == (lusifer_cdensity.nsinv[0][invariants.size()] - 1) );

                invariants.emplace_back(b, a);
            }

            assert( index == lusifer_cdensity.numinv[0][a][b] );

            channels_.at(a).invariants.at(b).index = index - 1;
        }
    }

    invariants.shrink_to_fit();

    counter = 0;

    for (std::size_t a = 0; a != channels_.size(); ++a)
    {
        for (std::size_t b = 0; b != channels_.at(a).processes.size(); ++b)
        {
            std::size_t index = 0;

            for (std::size_t c = 0; c != a; ++c)
            {
                for (std::size_t d = 0; d != channels_.at(c).processes.size(); ++d)
                {
                    if (processes_equal(a, b, c, d, nex))
                    {
                        index = channels_.at(c).processes.at(d).index + 1;

                        // break out of two loops
                        c = a - 1;
                        break;
                    }
                }
            }

            if (index == 0)
            {
                ++counter;
                index = counter;

                assert( a == (lusifer_cdensity.chprocess[0][processes_.size()] - 1) );
                assert( b == (lusifer_cdensity.ntprocess[0][processes_.size()] - 1) );

                processes_.emplace_back(b, a);
            }

            assert( index == lusifer_cdensity.numprocess[0][a][b] );

            channels_.at(a).processes.at(b).index = index - 1;
        }
    }

    processes_.shrink_to_fit();

    counter = 0;

    for (std::size_t a = 0; a != channels_.size(); ++a)
    {
        for (std::size_t b = 0; b != channels_.at(a).decays.size(); ++b)
        {
            std::size_t index = 0;

            for (std::size_t c = 0; c != a; ++c)
            {
                for (std::size_t d = 0; d != channels_.at(c).decays.size(); ++d)
                {
                    if (channels_.at(a).decays.at(b) == channels_.at(c).decays.at(d))
                    {
                        index = channels_.at(c).decays.at(d).index + 1;

                        // break out of two loops
                        c = a - 1;
                        break;
                    }
                }
            }

            if (index == 0)
            {
                ++counter;
                index = counter;

                assert( a == (lusifer_cdensity.chdecay[0][decays.size()] - 1) );
                assert( b == (lusifer_cdensity.nsdecay[0][decays.size()] - 1) );

                decays.emplace_back(b, a);
            }

            assert( index == lusifer_cdensity.numdecay[0][a][b] );

            channels_.at(a).decays.at(b).index = index - 1;
        }
    }

    decays.shrink_to_fit();

    // TODO: is the following needed?
    mcut.assign(std::begin(lusifer_cinv.mcutinv[0]), std::end(lusifer_cinv.mcutinv[0]));

    invariant_jacobians.reserve(invariants.size());
    process_jacobians.reserve(processes_.size());
    decay_jacobians.reserve(decays.size());
}

template <typename T>
std::size_t lusifer_psg<T>::channels() const
{
    return channels_.size();
}

template <typename T>
T lusifer_psg<T>::densities(std::vector<T>& densities)
{
    using std::acos;
    using std::fabs;
    using std::pow;
    using std::sqrt;

    std::size_t const allbinary = (1 << particles) - 1;

    // skip first three and last four elements; they are computed already
    for (std::size_t i = 0; i != allbinary - 1u; ++i)
    {
        // check if this isn't already computed
        if (p[i][0] != T())
        {
            continue;
        }

        // index of the sum that can be used with momentum conservation
        std::size_t const index = allbinary - i - 2;

        // if not already calculated, calculate all remaining sums and
        // invariants that are needed during the calculation of the densities
        // which have to be evaluated for all channels
        if (p[index][0] == T())
        {
            std::bitset<32> bitset_i(i + 1);

            std::size_t j = 0;
            while (!bitset_i.test(j))
            {
                ++j;
            }
            std::size_t leading_zero = 1 << j;

            p[i][0] = p[i - leading_zero][0] + p[leading_zero - 1][0];
            p[i][1] = p[i - leading_zero][1] + p[leading_zero - 1][1];
            p[i][2] = p[i - leading_zero][2] + p[leading_zero - 1][2];
            p[i][3] = p[i - leading_zero][3] + p[leading_zero - 1][3];

            // calculate the corresponding invariant
            s[i] = p[i][0] * p[i][0] - p[i][1] * p[i][1] - p[i][2] * p[i][2] - p[i][3] * p[i][3];

            p[index][0] = -p[i][0];
            p[index][1] = -p[i][1];
            p[index][2] = -p[i][2];
            p[index][3] = -p[i][3];

            s[index] = s[i];
        }
        else
        {
            // use momentum conservation
            p[i][0] = -p[index][0];
            p[i][1] = -p[index][1];
            p[i][2] = -p[index][2];
            p[i][3] = -p[index][3];

            s[i] = s[index];
        }
    }

    for (auto const& info : invariants)
    {
        auto const& invariant = channels_.at(info.channel).invariants.at(info.index);

        // index corresponding to the ns'th invariant
        std::size_t inv1 = invariant.in;
        // `inv1 + inv2 = allbinary - 3`
        std::size_t inv2 = (allbinary - 3 - inv1) - 2;

        T mmin = mcut.at(inv1);
        T mmax = mcut.at(inv2);

        for (std::size_t i = 0; &channels_[info.channel].invariants[i] != &invariant; ++i)
        {
            std::size_t const virt = channels_[info.channel].invariants[i].in;
            bool const condition = s[virt] > mcut[virt] * mcut[virt];

            // is there a minimum limit on this invariant?
            if (invariant.lmin.test(i) && condition)
            {
                // TODO: is `>` the right condition here?
                std::size_t const inv3 = inv1 - virt;
                mmin += sqrt(s[virt]) - mcut.at(inv1) + mcut.at(inv3);
                inv1 = inv3;
            }

            // is there a maximum limit on this invariant?
            if (invariant.lmax.test(i) && condition)
            {
                std::size_t const inv3 = inv2 - virt;
                mmax += sqrt(s[virt]) - mcut.at(inv2) + mcut.at(inv3);
                inv2 = inv3;
            }
        }

        T const smin = mmin * mmin;
        T const smax = (cmf_energy_ - mmax) * (cmf_energy_ - mmax);

        invariant_jacobians.push_back(jacobian(
            particle_infos.at(invariant.idhep).power,
            particle_infos.at(invariant.idhep).mass,
            particle_infos.at(invariant.idhep).width,
            s[invariant.in],
            smin,
            smax
        ));
    }

    for (auto const& info : processes_)
    {
        auto const& process = channels_[info.channel].processes[info.index];

        T const s  = this->s[process.in];
        T const s1 = this->s[process.out1];
        T const s2 = this->s[process.out2];
        T const t  = this->s[process.virt];
        T const t1 = this->s[process.in1];
        T const t2 = this->s[process.in2];

        auto const& tinv = calc_tinv(s, s1, s2, t1, t2);
        T const factor = T(2.0) * tinv.lambdat / acos(T(-1.0));

        process_jacobians.push_back(factor * jacobian(
             particle_infos.at(process.idhep).power,
            -particle_infos.at(process.idhep).mass,
             particle_infos.at(process.idhep).width,
            -t,
            -tinv.tmax,
            -tinv.tmin
        ));
    }

    for (auto const& info : decays)
    {
        auto const& decay = channels_[info.channel].decays[info.index];

        T const s  = this->s[decay.in];
        T const s1 = this->s[decay.out1];
        T const s2 = this->s[decay.out2];

        T const jacobian = T(2.0) * s / (acos(T(-1.0)) *
            hep::sqrt_kaellen(s, s1, s2));

        decay_jacobians.push_back(jacobian);
    }

    densities.clear();

    invariant_jacobians.clear();
    process_jacobians.clear();
    decay_jacobians.clear();

    // compute the jacobians for each channel
    for (auto const& channel : channels_)
    {
        T jacobian = T(1.0);

        for (auto const& invariant : channel.invariants)
        {
            jacobian *= invariant_jacobians[invariant.index];
        }

        for (auto const& process : channel.processes)
        {
            jacobian *= process_jacobians[process.index];
        }

        for (auto const& decay : channel.decays)
        {
            jacobian *= decay_jacobians[decay.index];
        }

        densities.push_back(jacobian);
    }

    return pow(T(0.5) / acos(T(-1.0)), T(3 * particles - 10));
}

template <typename T>
std::size_t lusifer_psg<T>::dimensions() const
{
    return 3 * (particles - 4) + 2 + extra_random_numbers;
}

template <typename T>
void lusifer_psg<T>::generate(
    std::vector<T> const& random_numbers,
    std::vector<T>& momenta,
    T cmf_energy,
    std::size_t channel
) {
    using std::acos;
    using std::atan2;
    using std::copysign;
    using std::fabs;
    using std::fmax;
    using std::fmin;
    using std::sqrt;

    for (auto& momentum : p)
    {
        momentum[0] = T();
    }

    auto r = random_numbers.begin();

    std::size_t const allbinary = (1 << particles) - 1;

    // first momentum: negative of the first beam momentum
    p[0] = { T(-0.5) * cmf_energy, T(), T(), T(0.5) * cmf_energy };
    // second momentum: negative of the second beam momentum
    p[1] = { T(-0.5) * cmf_energy, T(), T(), T(-0.5) * cmf_energy };
    // sum of the first two momenta
    p[2] = { T(-1.0) * cmf_energy, T(), T(), T() };

    // sum of all momenta except the first two
    p[allbinary - 4] = { cmf_energy, T(), T(), T() };
    // sum of all momenta except the second
    p[allbinary - 3] = { T( 0.5) * cmf_energy, T(), T(), T(0.5) * cmf_energy };
    // sum of all momenta except the first
    p[allbinary - 2] = { T( 0.5) * cmf_energy, T(), T(), T(-0.5) * cmf_energy };

    // square of the sum of the first two momenta
    s[2] = cmf_energy * cmf_energy;

    // square of the sum of all momenta except the first two
    s[allbinary - 4] = s[2];

    // iteratively construct the time-like invariants
    for (auto const& invariant : channels_[channel].invariants)
    {
        // index corresponding to the ns'th invariant
        std::size_t inv1 = invariant.in;
        // `inv1 + inv2 = allbinary - 3`
        std::size_t inv2 = (allbinary - 3 - inv1) - 2;

        T mmin = mcut[inv1];
        T mmax = mcut[inv2];

        for (std::size_t i = 0; &channels_[channel].invariants[i] != &invariant; ++i)
        {
            std::size_t const virt = channels_[channel].invariants[i].in;
            bool const condition = s[virt] > mcut[virt] * mcut[virt];

            // is there a minimum limit on this invariant?
            if (invariant.lmin.test(i) && condition)
            {
                std::size_t const inv3 = inv1 - virt;
                mmin += sqrt(s[virt]) - mcut[inv1] + mcut[inv3];
                inv1 = inv3;
            }

            // is there a maximum limit on this invariant?
            if (invariant.lmax.test(i) && condition)
            {
                std::size_t const inv3 = inv2 - virt;
                mmax += sqrt(s[virt]) - mcut[inv2] + mcut[inv3];
                inv2 = inv3;
            }
        }

        T const smin = mmin * mmin;
        T const smax = (cmf_energy - mmax) * (cmf_energy - mmax);

        // TODO: replace particle_infos calls with model

        s[invariant.in] = fabs(map(
            particle_infos.at(invariant.idhep).power,
            particle_infos.at(invariant.idhep).mass,
            particle_infos.at(invariant.idhep).width,
            *r++,
            smin,
            smax
        ));
    }

    for (auto const& process : channels_[channel].processes)
    {
        T const s  = this->s[process.in];
        T const s1 = this->s[process.out1];
        T const s2 = this->s[process.out2];
        T const t1 = this->s[process.in1];
        T const t2 = this->s[process.in2];
        T& t = this->s[process.virt];

        auto const& tinv = calc_tinv(s, s1, s2, t1, t2);
        T phi = T(2.0) * acos(T(-1.0)) * *r++;

        t = -map(
             particle_infos.at(process.idhep).power,
            -particle_infos.at(process.idhep).mass,
             particle_infos.at(process.idhep).width,
            *r++,
            -tinv.tmax,
            -tinv.tmin
        );

        T cos_theta = ((s + s1 - s2) * (s + t1 - t2) - T(2.0) * s *
            (t1 + s1 - t)) / (tinv.lambdas * tinv.lambdat);

        // make sure the cosine is within the proper bounds
        cos_theta = fmax(T(-1.0), fmin(cos_theta, T(1.0)));

        auto& p1 = p[process.out1];

        T const sqrts = sqrt(s);

        p1 = { T(0.5) * (s + s1 - s2) / sqrts, T(), T(), T(0.5) * tinv.lambdas / sqrts };

        auto const& q1 = p[process.in1];
        phi = copysign(phi, q1[3]);

        rotate(p1, phi, cos_theta);

        auto const& q2 = p[process.in2];

        std::array<T, 4> const q = { q1[0] + q2[0], q1[1] + q2[1], q1[2] + q2[2], q1[3] + q2[3] };

        std::array<T, 4> k1 = { q1[0], q1[1], q1[2], q1[3] };

        hep::boost(sqrts, k1, q, false);

        // the following can not be abbreviated with a single call to atan2, is
        // this choice of phi consistent?
        if (k1[1] == T())
        {
            phi = copysign(T(0.5) * acos(T(-1.0)), k1[2]);
        }
        else
        {
            phi = atan2(k1[2], k1[1]);
        }

        cos_theta = k1[3] / sqrt(k1[1] * k1[1] + k1[2] * k1[2] + k1[3] * k1[3]);

        auto& p2 = p[process.out2];
        decay_momenta(sqrts, q, -phi, cos_theta, p1, p2);

        auto& qt = p[process.virt];
        qt = { q1[0] - p1[0], q1[1] - p1[1], q1[2] - p1[2], q1[3] - p1[3] };
    }

    for (auto const& decay : channels_[channel].decays)
    {
        T const s  = this->s[decay.in];
        T const s1 = this->s[decay.out1];
        T const s2 = this->s[decay.out2];

        auto& p1 = p[decay.out1];
        T const sqrts = sqrt(s);

        p1 = { T(0.5) * (s + s1 - s2) / sqrts, T(), T(),
            T(0.5) * hep::sqrt_kaellen(s, s1, s2) / sqrts };

        T const phi = T(2.0) * acos(T(-1.0)) * *r++;
        T const cos_theta = T(2.0) * *r++ - T(1.0);

        decay_momenta(sqrts, p[decay.in], phi, cos_theta, p1, p[decay.out2]);
    }

    // assign beam momenta
    momenta.at(0) = T( 0.5) * cmf_energy;
    momenta.at(1) = T();
    momenta.at(2) = T();
    momenta.at(3) = T(-0.5) * cmf_energy;
    momenta.at(4) = T( 0.5) * cmf_energy;
    momenta.at(5) = T();
    momenta.at(6) = T();
    momenta.at(7) = T( 0.5) * cmf_energy;

    for (std::size_t i = 2; i != particles; ++i)
    {
        momenta.at(4 * i + 0) = p[(1 << i) - 1][0];
        momenta.at(4 * i + 1) = p[(1 << i) - 1][1];
        momenta.at(4 * i + 2) = p[(1 << i) - 1][2];
        momenta.at(4 * i + 3) = p[(1 << i) - 1][3];
    }

    cmf_energy_ = cmf_energy;

    // append the extra remaining random numbers to the end of `momenta`
    std::copy_n(r, extra_random_numbers,
        std::prev(momenta.end(), extra_random_numbers));
}

template <typename T>
std::size_t lusifer_psg<T>::map_dimensions() const
{
    return 4 * particles + extra_random_numbers;
}

}

namespace hep
{

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_lusifer_phase_space_generator(
    T min_energy,
    T cmf_energy,
    std::vector<std::string> const& processes,
    lusifer_constants<T> const& constants,
    std::size_t extra_random_numbers
) {
    return std::make_unique<hadron_hadron_psg_adapter<lusifer_psg<T>>>(
        min_energy,
        cmf_energy,
        processes,
        constants,
        extra_random_numbers
    );
}

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_lusifer_phase_space_generator(
    T min_energy,
    T cmf_energy,
    std::string const& process,
    lusifer_constants<T> const& constants,
    std::size_t extra_random_numbers
) {
    return make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        std::vector<std::string>{process},
        constants,
        extra_random_numbers
    );
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template std::unique_ptr<phase_space_generator<double>>
make_lusifer_phase_space_generator(
    double,
    double,
    std::vector<std::string> const&,
    lusifer_constants<double> const&,
    std::size_t
);

template std::unique_ptr<phase_space_generator<double>>
make_lusifer_phase_space_generator(
    double,
    double,
    std::string const&,
    lusifer_constants<double> const&,
    std::size_t
);

}
