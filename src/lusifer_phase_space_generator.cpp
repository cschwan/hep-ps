#include "hep/ps/fortran_helper.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "hep/ps/lusifer_ps_channels.hpp"
#include "hep/ps/lusifer_ps_functions.hpp"
#include "hep/ps/ps_functions.hpp"

#include "hadron_hadron_psg_adapter.hpp"
// FIXME: this header should not be needed
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
    std::bitset<16 * 8> upper_bound;
    std::bitset<16 * 8> lower_bound;
    std::size_t index;
};

struct inequality
{
    inequality(std::size_t lhs, std::size_t rhs1, std::size_t rhs2)
        : lhs(lhs)
        , rhs1(rhs1)
        , rhs2(rhs2)
    {
    }

    std::size_t lhs;
    std::size_t rhs1;
    std::size_t rhs2;
};

inequality get_inequality_for_invariant(
    std::vector<invariant> const& invariants,
    std::size_t inv_index
) {
    auto const in = inv_index + 1;

    // for invariants containing two momenta, the invariants are the external momenta squared
    if (std::bitset<32>(in).count() == 2)
    {
        auto const out1 = in & -in;
        auto const out2 = (in ^ out1) & -(in ^ out1);

        return { in - 1, out1 - 1, out2 - 1 };
    }

    assert( !invariants.empty() );

    for (auto i = invariants.begin(); i != invariants.end(); ++i)
    {
        auto const in1 = i->in + 1;

        if (std::bitset<32>(in ^ in1).count() == 1)
        {
            return { in - 1, in1 - 1, (in ^ in1) - 1 };
        }

        for (auto j = std::next(i); j != invariants.end(); ++j)
        {
            auto const in2 = j->in + 1;

            if (((in1 | in2) == in) && ((in1 & in2) == 0))
            {
                return { in - 1, in1 - 1, in2 - 1 };
            }
        }
    }

    assert( false );
}

std::vector<inequality> build_inequalities(
    std::vector<invariant> const& invariants,
    std::size_t cmf_inv
) {
    std::vector<inequality> result;
    result.reserve(invariants.size() + 1);

    // first inequality from the center-of-mass energy
    result.push_back(get_inequality_for_invariant(invariants, cmf_inv));

    for (auto const& inv : invariants)
    {
        result.push_back(get_inequality_for_invariant(invariants, inv.in));
    }

    return result;
}

template <typename T>
std::bitset<8 * 16> integral_bound(
    std::vector<inequality>::const_iterator inv,
    std::vector<T> const& mcut,
    std::vector<std::size_t>& stack,
    std::vector<inequality> const& inequalities
) {
    std::vector<std::size_t> invs;
    invs.reserve(inequalities.size());

    while (!stack.empty())
    {
        std::size_t const index = stack.back();

        // remove invariant
        stack.pop_back();

        // if there is one single bit, it's an external mass
        if (std::bitset<64>(index + 1).count() == 1)
        {
            // only add the invariant if it's nonzero
            if (mcut.at(index + 1) != T())
            {
                invs.push_back(index);
            }
        }
        else // if it's not an external mass, it is an invariant
        {
            auto const result = std::find_if(inequalities.begin(), inequalities.end(),
                [=](inequality const& ineq) { return ineq.lhs == index; });

            assert( result != inequalities.end() );

            // is the invariant not generated yet?
            if (std::distance(inv, result) >= 0)
            {
                stack.push_back(result->rhs1);
                stack.push_back(result->rhs2);
            }
            else
            {
                invs.push_back(index);
            }
        }
    }

    // sort invariants to optimize access pattern in later calculations
    std::sort(invs.begin(), invs.end());
    std::bitset<8 * 16> bound;

    for (std::size_t i = 0; i != invs.size(); ++i)
    {
        bound |= std::bitset<8 * 16>(invs.at(i)) << (16 * i);
    }

    return bound;
}

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
void decay_momenta(
    T m,
    std::array<T, 4> const& q,
    T phi,
    T cos_theta,
    std::array<T, 4>& p1,
    std::array<T, 4>& p2
) {
    hep::rotate(p1, phi, cos_theta);
    hep::boost(m, p1, q, true);

    p2[0] = q[0] - p1[0];
    p2[1] = q[1] - p1[1];
    p2[2] = q[2] - p1[2];
    p2[3] = q[3] - p1[3];
}

template <typename T>
class lusifer_psg
{
public:
    using numeric_type = T;

    lusifer_psg(
        std::vector<hep::ps_channel> const& channels,
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

    if (a.lower_bound != b.lower_bound)
    {
        return false;
    }

    if (a.upper_bound != b.upper_bound)
    {
        return false;
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
    std::vector<hep::ps_channel> const& channels,
    hep::lusifer_constants<T> const& constants,
    std::size_t extra_random_numbers
) {
    auto const nex = channels.front().invariants().size() + 4;

    s.resize(1 << nex);
    p.resize(1 << nex);

    this->particles = nex;
    this->extra_random_numbers = extra_random_numbers;
    this->channels_.resize(channels.size());

    // FIXME: access to FORTRAN common block
    mcut.assign(std::begin(lusifer_cinv.mcutinv[0]), std::end(lusifer_cinv.mcutinv[0]));

    // TODO: assign external masses to `s` vector
    for (std::size_t i = 0; i != particles; ++i)
    {
        // FIXME: access to FORTRAN common block
        T const m = lusifer_cinv.mcutinv[0][1 << i];
        s[(1 << i) - 1] = m * m;
    }

    std::vector<std::size_t> inv_lo;
    std::vector<std::size_t> inv_hi;
    std::vector<std::size_t> stack;

    for (std::size_t i = 0; i != channels.size(); ++i)
    {
        auto& channel = this->channels_.at(i);

        channel.invariants.reserve(channels.at(i).invariants().size());
        channel.processes.reserve(channels.at(i).tchannels().size());
        channel.decays.reserve(channels.at(i).decays().size());

        for (std::size_t j = 0; j != channel.invariants.capacity(); ++j)
        {
            channel.invariants.emplace_back(
                channels.at(i).invariants().at(j).mom_id(),
                channels.at(i).invariants().at(j).pdg_id()
            );
        }

        std::size_t const allbinary = (1 << nex) - 1;

        auto const& inequalities = build_inequalities(channel.invariants, allbinary - 4);

        // skip first inequality, because it is already generated elsewhere
        for (auto i = std::next(inequalities.cbegin()); i != inequalities.cend(); ++i)
        {
            stack.clear();
            stack.push_back(i->lhs);

            std::size_t const index = std::distance(inequalities.cbegin(), i) - 1;
            channel.invariants.at(index).lower_bound = integral_bound(i, mcut, stack, inequalities);
        }

        // skip first inequality, because it is already generated elsewhere
        for (auto i = std::next(inequalities.cbegin()); i != inequalities.cend(); ++i)
        {
            stack.clear();

            std::size_t positive = i->lhs;
            std::vector<inequality>::const_iterator next;

            do
            {
                next = std::find_if(inequalities.begin(), inequalities.end(),
                    [=](inequality const& ineq) {
                        return (ineq.rhs1 == positive) || (ineq.rhs2 == positive);
                });

                stack.push_back((positive == next->rhs1) ? next->rhs2 : next->rhs1);
                positive = next->lhs;
            }
            while (std::distance(i, next) >= 0);

            auto upper_bound_neg = integral_bound(i, mcut, stack, inequalities);
            std::size_t const index = std::distance(inequalities.cbegin(), i) - 1;

            channel.invariants.at(index).upper_bound =
                (upper_bound_neg << 16) | std::bitset<8 * 16>(positive);
        }

        for (std::size_t j = 0; j != channel.processes.capacity(); ++j)
        {
            channel.processes.emplace_back(
                channels.at(i).tchannels().at(j).in1_mom_id(),
                channels.at(i).tchannels().at(j).in2_mom_id(),
                channels.at(i).tchannels().at(j).out1_mom_id(),
                channels.at(i).tchannels().at(j).out2_mom_id(),
                channels.at(i).tchannels().at(j).in_mom_id(),
                channels.at(i).tchannels().at(j).virt_id(),
                channels.at(i).tchannels().at(j).pdg_id()
            );
        }

        for (std::size_t j = 0; j != channel.decays.capacity(); ++j)
        {
            channel.decays.emplace_back(
                channels.at(i).decays().at(j).in_mom_id(),
                channels.at(i).decays().at(j).out1_mom_id(),
                channels.at(i).decays().at(j).out2_mom_id()
            );
        }
    }

    particle_infos.reserve(maxv + 1);

    for (std::size_t i = 0; i != particle_infos.capacity(); ++i)
    {
        // TODO: make this parameter available from outside
        T power = T(0.9);
        // FIXME: access to FORTRAN common block
        T mass  = lusifer_general.mass[i];
        // FIXME: access to FORTRAN common block
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

                invariants.emplace_back(b, a);
            }

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

                processes_.emplace_back(b, a);
            }

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

                decays.emplace_back(b, a);
            }

            channels_.at(a).decays.at(b).index = index - 1;
        }
    }

    decays.shrink_to_fit();


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

        T min = T();

        auto lower = invariant.lower_bound;
        auto upper = invariant.upper_bound;

        while (lower.any())
        {
            auto const index = (lower & std::bitset<8 * 16>((1 << 16) - 1)).to_ulong();
            lower >>= 16;
            min += sqrt(s.at(index));
        }

        assert( upper.any() );

        auto const index = (upper & std::bitset<8 * 16>((1 << 16) - 1)).to_ulong();
        upper >>= 16;
        T max_pos = (index == allbinary - 4) ? cmf_energy_ : sqrt(s.at(index));
        T max_neg = T();

        while (upper.any())
        {
            auto const index = (upper & std::bitset<8 * 16>((1 << 16) - 1)).to_ulong();
            upper >>= 16;
            max_neg += sqrt(s.at(index));
        }

        T const max = max_pos - max_neg;
        T const smin = min * min;
        T const smax = max * max;

        invariant_jacobians.push_back(hep::lusifer_invariant_jacobian(
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

        auto const& tinv = hep::calculate_space_like_invariant_bounds(s, s1, s2, t1, t2);
        T const factor = T(2.0) * tinv.lambdat / acos(T(-1.0));

        process_jacobians.push_back(factor * hep::lusifer_invariant_jacobian(
             particle_infos.at(process.idhep).power,
            -particle_infos.at(process.idhep).mass,
             particle_infos.at(process.idhep).width,
            -t,
            -tinv.max,
            -tinv.min
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
        T min = T();

        auto lower = invariant.lower_bound;
        auto upper = invariant.upper_bound;

        while (lower.any())
        {
            auto const index = (lower & std::bitset<8 * 16>((1 << 16) - 1)).to_ulong();
            lower >>= 16;
            min += sqrt(s.at(index));
        }

        assert( upper.any() );

        auto const index = (upper & std::bitset<8 * 16>((1 << 16) - 1)).to_ulong();
        upper >>= 16;
        T max_pos = (index == allbinary - 4) ? cmf_energy : sqrt(s.at(index));
        T max_neg = T();

        while (upper.any())
        {
            auto const index = (upper & std::bitset<8 * 16>((1 << 16) - 1)).to_ulong();
            upper >>= 16;
            max_neg += sqrt(s.at(index));
        }

        T const max = max_pos - max_neg;
        T const smin = min * min;
        T const smax = max * max;

        // TODO: replace particle_infos calls with model

        s[invariant.in] = fabs(hep::lusifer_invariant_map(
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

        auto const& tinv = hep::calculate_space_like_invariant_bounds(s, s1, s2, t1, t2);
        T phi = T(2.0) * acos(T(-1.0)) * *r++;

        t = -hep::lusifer_invariant_map(
             particle_infos.at(process.idhep).power,
            -particle_infos.at(process.idhep).mass,
             particle_infos.at(process.idhep).width,
            *r++,
            -tinv.max,
            -tinv.min
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

        hep::rotate(p1, phi, cos_theta);

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
    std::vector<ps_channel> const& channels,
    lusifer_constants<T> const& constants,
    std::size_t extra_random_numbers
) {
    return std::make_unique<hadron_hadron_psg_adapter<lusifer_psg<T>>>(
        min_energy,
        cmf_energy,
        channels,
        constants,
        extra_random_numbers
    );
}

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_lusifer_phase_space_generator(
    T min_energy,
    T cmf_energy,
    std::vector<std::string> const& processes,
    lusifer_constants<T> const& constants,
    std::size_t extra_random_numbers
) {
    return make_lusifer_phase_space_generator(
        min_energy,
        cmf_energy,
        lusifer_ps_channels(processes, constants),
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
    std::vector<ps_channel> const&,
    lusifer_constants<double> const&,
    std::size_t
);

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
