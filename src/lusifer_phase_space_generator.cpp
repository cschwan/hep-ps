#include "hep/ps/boost.hpp"
#include "hep/ps/fortran_helper.hpp"
#include "hep/ps/kaellen.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"

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

struct invariant
{
    invariant(
        std::size_t in,
        std::size_t idhep,
        std::bitset<std::numeric_limits<std::size_t>::digits> lmin,
        std::bitset<std::numeric_limits<std::size_t>::digits> lmax,
        std::size_t index
    )
        : in(in)
        , idhep(idhep)
        , lmin(lmin)
        , lmax(lmax)
        , index(index)
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
        std::size_t idhep,
        std::size_t index
    )
        : in1(in1)
        , in2(in2)
        , out1(out1)
        , out2(out2)
        , in(in)
        , virt(virt)
        , idhep(idhep)
        , index(index)
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
        std::size_t out2,
        std::size_t index
    )
        : in(in)
        , out1(out1)
        , out2(out2)
        , index(index)
    {
    }

    std::size_t in;
    std::size_t out1;
    std::size_t out2;
    std::size_t index;
};

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
    std::size_t max_particles;
    std::size_t extra_random_numbers;
    T cmf_energy_;
};

template <typename T>
lusifer_psg<T>::lusifer_psg(
    std::vector<std::string> const& processes,
    hep::lusifer_constants<T> const& constants,
    std::size_t extra_random_numbers
) {
    int maxex;
    int maxgen;
    lusifer_extra_max(&maxex, &maxgen);

    int nex = 0;

    // FORTRAN counting: use the first generator
    int g = 1;

    for (auto const& process : processes)
    {
        // each particle must be specified with three characters
        assert( process.size() % 3 == 0 );

        if (nex == 0)
        {
            // there must be at least four particles
            assert( process.size() >= 4 * 3 );

            nex = process.size() / 3;

            double mw = static_cast <double> (constants.mass_w);
            double gw = static_cast <double> (constants.width_w);
            double mz = static_cast <double> (constants.mass_z);
            double gz = static_cast <double> (constants.width_z);
            double mh = static_cast <double> (constants.mass_h);
            double gh = static_cast <double> (constants.width_h);
            double mt = static_cast <double> (constants.mass_t);
            double gt = static_cast <double> (constants.width_t);

            // set constants and the number of particles
            lusifer_extra_set(g, nex, mw, gw, mz, gz, mh, gh, mt, gt);

            // the number of particles must be supported by the generator
            assert( nex <= maxex );
        }

        // all processes must have the same number of particles
        assert( nex == static_cast <int> (process.size() / 3) );

        // fill up the string with three spaces for particle not used
        std::string process0 = process;
        process0.append(3 * (maxex - nex), ' ');

        // TODO: what is the meaning of this parameter?
        int lightfermions = 0;
        // do not include cuts in the phase space generation
        int includecuts = 0;
        // do not print channel information
        int sout = 0;

        lusifer_initphasespace(
            process0.c_str(),
            g,
            lightfermions,
            includecuts,
            sout
        );
    }

    int channels;
    lusifer_extra_data(g, &channels);

    // there must be at least one channel, otherwise something went wrong
    assert( channels > 0 );

    s.resize(1 << nex);
    p.resize(1 << nex);

    this->particles = nex;
    this->max_particles = maxe;
    this->extra_random_numbers = extra_random_numbers;
    this->channels_.resize(channels);

    for (std::size_t i = 0; i != static_cast <std::size_t> (channels); ++i)
    {
        auto& channel = this->channels_.at(i);

        channel.invariants.reserve(lusifer_cinv.ninv[0][i]);
        channel.processes.reserve(lusifer_cprocess.nprocess[0][i]);
        channel.decays.reserve(lusifer_cdecay.ndecay[0][i]);

        for (std::size_t j = 0; j != channel.invariants.capacity(); ++j)
        {
            decltype (invariant::lmin) lmin;
            decltype (invariant::lmax) lmax;

            for (std::size_t k = 0; k != maxe; ++k)
            {
                lmin.set(k, lusifer_cinv.lmin[0][i][j][k]);
                lmax.set(k, lusifer_cinv.lmax[0][i][j][k]);
            }

            channel.invariants.emplace_back(
                lusifer_cinv.ininv[0][i][j] - 1,
                lusifer_cinv.idhepinv[0][i][j],
                lmin,
                lmax,
                lusifer_cdensity.numinv[0][i][j] - 1
            );
        }

        for (std::size_t j = 0; j != channel.processes.capacity(); ++j)
        {
            channel.processes.emplace_back(
                lusifer_cprocess.in1process[0][i][j] - 1,
                lusifer_cprocess.in2process[0][i][j] - 1,
                lusifer_cprocess.out1process[0][i][j] - 1,
                lusifer_cprocess.out2process[0][i][j] - 1,
                lusifer_cprocess.inprocess[0][i][j] - 1,
                lusifer_cprocess.virtprocess[0][i][j] - 1,
                lusifer_cprocess.idhepprocess[0][i][j],
                lusifer_cdensity.numprocess[0][i][j] - 1
            );
        }

        for (std::size_t j = 0; j != channel.decays.capacity(); ++j)
        {
            channel.decays.emplace_back(
                lusifer_cdecay.indecay[0][i][j] - 1,
                lusifer_cdecay.out1decay[0][i][j] - 1,
                lusifer_cdecay.out2decay[0][i][j] - 1,
                lusifer_cdensity.numdecay[0][i][j] - 1
            );
        }
    }

    invariants.reserve(lusifer_cdensity.maxinv[0]);

    for (std::size_t i = 0; i != invariants.capacity(); ++i)
    {
        invariants.emplace_back(
            lusifer_cdensity.nsinv[0][i] - 1,
            lusifer_cdensity.chinv[0][i] - 1
        );
    }

    processes_.reserve(lusifer_cdensity.maxprocess[0]);

    for (std::size_t i = 0; i != processes_.capacity(); ++i)
    {
        processes_.emplace_back(
            lusifer_cdensity.ntprocess[0][i] - 1,
            lusifer_cdensity.chprocess[0][i] - 1
        );
    }

    decays.reserve(lusifer_cdensity.maxdecay[0]);

    for (std::size_t i = 0; i != decays.capacity(); ++i)
    {
        decays.emplace_back(
            lusifer_cdensity.nsdecay[0][i] - 1,
            lusifer_cdensity.chdecay[0][i] - 1
        );
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
