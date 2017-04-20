#include "hep/ps/fortran_helper.hpp"
#include "hep/ps/kaellen.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "lusifer_interfaces.hpp"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iterator>
#include <limits>
#include <utility>

constexpr std::size_t maxg = 1;
constexpr std::size_t maxe = 9;
constexpr std::size_t maxch = 20000;
constexpr std::size_t maxv = 40;

extern "C"
{

struct
{
	double alphaisr;
	double scale;
	double meisr;
	double s[1 << maxe];
	double p[1 << maxe][4];
	double mass[maxv + 1];
	double width[maxv + 1];
	int nchannel[maxg];
	int nexternal[maxg];
	int allbinary[maxg];
} lusifer_general_;

struct
{
	double powerinv[maxg][maxch][maxe];
	double mcutinv[maxg][(1 << maxe) + 1];
	int ininv[maxg][maxch][maxe];
	int idhepinv[maxg][maxch][maxv];
	int ninv[maxg][maxch];
	int lmin[maxg][maxch][maxe][maxe];
	int lmax[maxg][maxch][maxe][maxe];
} lusifer_cinv_;

struct
{
	int indecay[maxg][maxch][maxe];
	int out1decay[maxg][maxch][maxe];
	int out2decay[maxg][maxch][maxe];
	int ndecay[maxg][maxch];
} lusifer_cdecay_;

struct
{
	double powerprocess[maxg][maxch][maxe];
	double ccutprocess[maxg][1 << maxe];
	int in1process[maxg][maxch][maxe];
	int in2process[maxg][maxch][maxe];
	int out1process[maxg][maxch][maxe];
	int out2process[maxg][maxch][maxe];
	int inprocess[maxg][maxch][maxe];
	int virtprocess[maxg][maxch][maxe];
	int idhepprocess[maxg][maxch][maxe];
	int nprocess[maxg][maxch];
} lusifer_cprocess_;

struct
{
	int nsinv[maxg][maxch];
	int chinv[maxg][maxch];
	int maxinv[maxg];
	int ntprocess[maxg][maxch];
	int chprocess[maxg][maxch];
	int maxprocess[maxg];
	int nsdecay[maxg][maxch];
	int chdecay[maxg][maxch];
	int maxdecay[maxg];
	int numinv[maxg][maxch][maxe];
	int numprocess[maxg][maxch][maxe];
	int numdecay[maxg][maxch][maxe];
} lusifer_cdensity_;

}

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
		T const min = atan((xmin - m2) / mg);
		T const max = atan((xmax - m2) / mg);
		T const xprime = x - m2;

		return mg / ((max - min) * (xprime * xprime + mg * mg));
	}

	if (power == T())
	{
		return T(1.0) / (xmax - xmin);
	}

	// TODO: make this parameter accessible from outside
	m2 -= T(1e-6);

	if (power == T(1.0))
	{
		// TODO: WARNING this branch is untested!
		return T(1.0) / (log((xmax - m2) / (xmin - m2)) * (x - m2));
	}

	T const omp = T(1.0) - power;

	return omp / ((pow(xmax - m2, omp) - pow(xmin - m2, omp)) *
		pow(x - m2, power));
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

	if (width > T())
	{
		T const mg = mass * width;
		T const min = atan((xmin - m2) / mg);
		T const max = atan((xmax - m2) / mg);

		return m2 + mg * tan(x * (max - min) + min);
	}

	if (power == T())
	{
		return x * xmax + (T(1.0) - x) * xmin;
	}

	// TODO: make this parameter accessible from outside
	m2 -= T(1e-6);

	if (power == T(1.0))
	{
		// TODO: WARNING this branch is untested!
		return m2 + exp(x * log(xmax - m2) + (T(1.0) - x) * log(xmin - m2));
	}

	T const omp = T(1.0) - power;

	return m2 + pow(x * pow(xmax - m2, omp) + (T(1.0) - x) *
		pow(xmin - m2, omp), T(1.0) / omp);
}

template <typename T>
void rotate(std::array<T, 4>& p, T phi, T cos_theta)
{
	using std::cos;
	using std::sin;
	using std::sqrt;

	T const sin_phi = sin(phi);
	T const cos_phi = cos(phi);
	T const sin_theta = sqrt((T(1.0) - cos_theta) * (T(1.0) + cos_theta));

	T const px = p[1];
	T const py = p[2];
	T const pz = p[3];

	p[1] =  px * cos_theta * cos_phi + py * sin_phi + pz * sin_theta * cos_phi;
	p[2] = -px * cos_theta * sin_phi + py * cos_phi - pz * sin_theta * sin_phi;
	p[3] = -px * sin_theta                          + pz * cos_theta;
}

template <typename T>
void boost(T m, std::array<T, 4>& p, std::array<T, 4> const& q, bool inverse)
{
	T const sign = inverse ? T(1.0) : T(-1.0);
	T const bx = sign * q[1] / m;
	T const by = sign * q[2] / m;
	T const bz = sign * q[3] / m;
	T const gamma = q[0] / m;
	T const a = T(1.0) / (T(1.0) + gamma);
	T const p0 = p[0];
	T const px = p[1];
	T const py = p[2];
	T const pz = p[3];
	T const bp = bx * px + by * py + bz * pz;

	p[0] = gamma * p0 + bp;
	p[1] = px + bx * p0 + a * bp * bx;
	p[2] = py + by * p0 + a * bp * by;
	p[3] = pz + bz * p0 + a * bp * bz;
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
	boost(m, p1, q, true);

	p2[0] = q[0] - p1[0];
	p2[1] = q[1] - p1[1];
	p2[2] = q[2] - p1[2];
	p2[3] = q[3] - p1[3];
}

}

namespace hep
{

template <typename T>
lusifer_constants<T>::lusifer_constants(
	T mass_h,
	T width_h,
	T mass_t,
	T width_t,
	T mass_w,
	T width_w,
	T mass_z,
	T width_z
)
	: mass_h(mass_h)
	, width_h(width_h)
	, mass_t(mass_t)
	, width_t(width_t)
	, mass_w(mass_w)
	, width_w(width_w)
	, mass_z(mass_z)
	, width_z(width_z)
{
}

template <typename T>
class lusifer_phase_space_generator<T>::impl
{
public:
	impl(
		std::size_t particles,
		std::size_t max_particles,
		std::size_t channels,
		std::size_t extra_random_numbers
	);

	std::vector<channel> channels;
	std::vector<invariant_info> invariants;
	std::vector<process_info> processes;
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
lusifer_phase_space_generator<T>::impl::impl(
	std::size_t particles,
	std::size_t max_particles,
	std::size_t channels,
	std::size_t extra_random_numbers
)
	: s(1 << particles)
	, p(1 << particles)
	, particles(particles)
	, max_particles(max_particles)
	, extra_random_numbers(extra_random_numbers)
{
	this->channels.resize(channels);

	for (std::size_t i = 0; i != channels; ++i)
	{
		auto& channel = this->channels.at(i);

		channel.invariants.reserve(lusifer_cinv_.ninv[0][i]);
		channel.processes.reserve(lusifer_cprocess_.nprocess[0][i]);
		channel.decays.reserve(lusifer_cdecay_.ndecay[0][i]);

		for (std::size_t j = 0; j != channel.invariants.capacity(); ++j)
		{
			decltype (invariant::lmin) lmin;
			decltype (invariant::lmax) lmax;

			for (std::size_t k = 0; k != maxe; ++k)
			{
				lmin.set(k, lusifer_cinv_.lmin[0][i][j][k]);
				lmax.set(k, lusifer_cinv_.lmax[0][i][j][k]);
			}

			channel.invariants.emplace_back(
				lusifer_cinv_.ininv[0][i][j] - 1,
				lusifer_cinv_.idhepinv[0][i][j],
				lmin,
				lmax,
				lusifer_cdensity_.numinv[0][i][j] - 1
			);
		}

		for (std::size_t j = 0; j != channel.processes.capacity(); ++j)
		{
			channel.processes.emplace_back(
				lusifer_cprocess_.in1process[0][i][j] - 1,
				lusifer_cprocess_.in2process[0][i][j] - 1,
				lusifer_cprocess_.out1process[0][i][j] - 1,
				lusifer_cprocess_.out2process[0][i][j] - 1,
				lusifer_cprocess_.inprocess[0][i][j] - 1,
				lusifer_cprocess_.virtprocess[0][i][j] - 1,
				lusifer_cprocess_.idhepprocess[0][i][j],
				lusifer_cdensity_.numprocess[0][i][j] - 1
			);
		}

		for (std::size_t j = 0; j != channel.decays.capacity(); ++j)
		{
			channel.decays.emplace_back(
				lusifer_cdecay_.indecay[0][i][j] - 1,
				lusifer_cdecay_.out1decay[0][i][j] - 1,
				lusifer_cdecay_.out2decay[0][i][j] - 1,
				lusifer_cdensity_.numdecay[0][i][j] - 1
			);
		}
	}

	invariants.reserve(lusifer_cdensity_.maxinv[0]);

	for (std::size_t i = 0; i != invariants.capacity(); ++i)
	{
		invariants.emplace_back(
			lusifer_cdensity_.nsinv[0][i] - 1,
			lusifer_cdensity_.chinv[0][i] - 1
		);
	}

	processes.reserve(lusifer_cdensity_.maxprocess[0]);

	for (std::size_t i = 0; i != processes.capacity(); ++i)
	{
		processes.emplace_back(
			lusifer_cdensity_.ntprocess[0][i] - 1,
			lusifer_cdensity_.chprocess[0][i] - 1
		);
	}

	decays.reserve(lusifer_cdensity_.maxdecay[0]);

	for (std::size_t i = 0; i != decays.capacity(); ++i)
	{
		decays.emplace_back(
			lusifer_cdensity_.nsdecay[0][i] - 1,
			lusifer_cdensity_.chdecay[0][i] - 1
		);
	}

	particle_infos.reserve(maxv + 1);

	for (std::size_t i = 0; i != particle_infos.capacity(); ++i)
	{
		// TODO: make this parameter available from outside
		T power = T(0.9);
		T mass  = lusifer_general_.mass[i];
		T width = lusifer_general_.width[i];

		if ((width > T()) || (i == 0) || ((i >= 30) && (i < 33)))
		{
			power = T();
		}

		particle_infos.emplace_back(power, mass, width);
	}

	// TODO: is the following needed?
	mcut.assign(
		std::begin(lusifer_cinv_.mcutinv[0]),
		std::end(lusifer_cinv_.mcutinv[0])
	);

	invariant_jacobians.reserve(invariants.size());
	process_jacobians.reserve(processes.size());
	decay_jacobians.reserve(decays.size());
}

template <typename T>
lusifer_phase_space_generator<T>::lusifer_phase_space_generator(
	std::vector<std::string> const& processes,
	lusifer_constants<T> const& constants,
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

			double mw = constants.mass_w;
			double gw = constants.width_w;
			double mz = constants.mass_z;
			double gz = constants.width_z;
			double mh = constants.mass_h;
			double gh = constants.width_h;
			double mt = constants.mass_t;
			double gt = constants.width_t;

			// set constants and the number of particles
			lusifer_extra_set(&g, &nex, &mw, &gw, &mz, &gz, &mh, &gh, &mt, &gt);

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
			&g,
			&lightfermions,
			&includecuts,
			&sout,
			process0.size()
		);
	}

	int channels;
	lusifer_extra_data(&g, &channels);

	// there must be at least one channel, otherwise something went wrong
	assert( channels > 0 );

	pimpl = std::move(std::unique_ptr<impl>(new impl(
		nex,
		maxe,
		channels,
		extra_random_numbers
	)));
}

template <typename T>
lusifer_phase_space_generator<T>::lusifer_phase_space_generator(
	std::string const& process,
	lusifer_constants<T> const& constants,
	std::size_t extra_random_numbers
)
	: lusifer_phase_space_generator(
		std::vector<std::string>{process},
		constants,
		extra_random_numbers
	)
{
}

template <typename T>
lusifer_phase_space_generator<T>::lusifer_phase_space_generator(
	lusifer_phase_space_generator<T>&& psg
)
	: pimpl(std::move(psg.pimpl))
{
}

template <typename T>
lusifer_phase_space_generator<T>::~lusifer_phase_space_generator() = default;

template <typename T>
std::size_t lusifer_phase_space_generator<T>::channels() const
{
	return pimpl->channels.size();
}

template <typename T>
T lusifer_phase_space_generator<T>::densities(std::vector<T>& densities)
{
	using std::acos;
	using std::fabs;
	using std::pow;
	using std::sqrt;

	auto& p = pimpl->p;
	auto& s = pimpl->s;

	std::size_t const allbinary = (1 << pimpl->particles) - 1;

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
			s[i] = p[i][0] * p[i][0] - p[i][1] * p[i][1] -
				p[i][2] * p[i][2] - p[i][3] * p[i][3];

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

	for (auto const& info : pimpl->invariants)
	{
		auto const& invariant = pimpl->channels.at(info.channel).invariants.at(info.index);

		// index corresponding to the ns'th invariant
		std::size_t inv1 = invariant.in;
		// `inv1 + inv2 = allbinary - 3`
		std::size_t inv2 = (allbinary - 3 - inv1) - 2;

		T mmin = pimpl->mcut.at(inv1);
		T mmax = pimpl->mcut.at(inv2);

		for (std::size_t i = 0;
			&pimpl->channels[info.channel].invariants[i] != &invariant; ++i)
		{
			std::size_t const virt = pimpl->channels[info.channel].invariants[i].in;
			bool const condition = s[virt] > pimpl->mcut[virt] * pimpl->mcut[virt];

			// is there a minimum limit on this invariant?
			if (invariant.lmin.test(i) && condition)
			{
				// TODO: is `>` the right condition here?
				std::size_t const inv3 = inv1 - virt;
				mmin += sqrt(s[virt]) - pimpl->mcut.at(inv1) + pimpl->mcut.at(inv3);
				inv1 = inv3;
			}

			// is there a maximum limit on this invariant?
			if (invariant.lmax.test(i) && condition)
			{
				std::size_t const inv3 = inv2 - virt;
				mmax += sqrt(s[virt]) - pimpl->mcut.at(inv2) + pimpl->mcut.at(inv3);
				inv2 = inv3;
			}
		}

		T const smin = mmin * mmin;
		T const smax = (pimpl->cmf_energy_ - mmax) * (pimpl->cmf_energy_ - mmax);

		pimpl->invariant_jacobians.push_back(jacobian(
			pimpl->particle_infos.at(invariant.idhep).power,
			pimpl->particle_infos.at(invariant.idhep).mass,
			pimpl->particle_infos.at(invariant.idhep).width,
			s[invariant.in],
			smin,
			smax
		));
	}

	for (auto const& info : pimpl->processes)
	{
		auto const& process = pimpl->channels[info.channel].processes[info.index];

		T const s  = pimpl->s[process.in];
		T const s1 = pimpl->s[process.out1];
		T const s2 = pimpl->s[process.out2];
		T const t  = pimpl->s[process.virt];
		T const t1 = pimpl->s[process.in1];
		T const t2 = pimpl->s[process.in2];

		T const lambdas = sqrt(kaellen(s, s1, s2));
		T const lambdat = sqrt(kaellen(s, t1, t2));

		T const tmp = (s + s1 - s2) * (s + t1 - t2);
		T const tmin = s1 + t1 - T(0.5) * (tmp + lambdas * lambdat) / s;
		T       tmax = s1 + t1 - T(0.5) * (tmp - lambdas * lambdat) / s;

		// TODO: make this parameter available from outside
		if (fabs(tmax) < T(1e-7))
		{
			tmax = T();
		}

		T const factor = T(2.0) * lambdat / acos(T(-1.0));

		pimpl->process_jacobians.push_back(factor * jacobian(
			 pimpl->particle_infos.at(process.idhep).power,
			-pimpl->particle_infos.at(process.idhep).mass,
			 pimpl->particle_infos.at(process.idhep).width,
			-t,
			-tmax,
			-tmin
		));
	}

	for (auto const& info : pimpl->decays)
	{
		auto const& decay = pimpl->channels[info.channel].decays[info.index];

		T const s  = pimpl->s[decay.in];
		T const s1 = pimpl->s[decay.out1];
		T const s2 = pimpl->s[decay.out2];

		T const jacobian = T(2.0) * s / (acos(T(-1.0)) *
			sqrt(kaellen(s, s1, s2)));

		pimpl->decay_jacobians.push_back(jacobian);
	}

	densities.clear();

	pimpl->invariant_jacobians.clear();
	pimpl->process_jacobians.clear();
	pimpl->decay_jacobians.clear();

	// compute the jacobians for each channel
	for (auto const& channel : pimpl->channels)
	{
		T jacobian = T(1.0);

		for (auto const& invariant : channel.invariants)
		{
			jacobian *= pimpl->invariant_jacobians[invariant.index];
		}

		for (auto const& process : channel.processes)
		{
			jacobian *= pimpl->process_jacobians[process.index];
		}

		for (auto const& decay : channel.decays)
		{
			jacobian *= pimpl->decay_jacobians[decay.index];
		}

		densities.push_back(jacobian);
	}

	return pow(T(0.5) / acos(T(-1.0)), T(3 * pimpl->particles - 10));
}

template <typename T>
std::size_t lusifer_phase_space_generator<T>::dimensions() const
{
	return 3 * (pimpl->particles - 4) + 2 + pimpl->extra_random_numbers;
}

template <typename T>
void lusifer_phase_space_generator<T>::generate(
	std::vector<T> const& random_numbers,
	std::vector<T>& momenta,
	T cmf_energy,
	std::size_t channel
) {
	using std::acos;
	using std::atan2;
	using std::copysign;
	using std::fabs;
	using std::sqrt;

	auto& s = pimpl->s;
	auto& p = pimpl->p;

	for (auto& momentum : p)
	{
		momentum[0] = T();
	}

	auto r = random_numbers.begin();

	std::size_t const allbinary = (1 << pimpl->particles) - 1;

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
	for (auto const& invariant : pimpl->channels[channel].invariants)
	{
		// index corresponding to the ns'th invariant
		std::size_t inv1 = invariant.in;
		// `inv1 + inv2 = allbinary - 3`
		std::size_t inv2 = (allbinary - 3 - inv1) - 2;

		T mmin = pimpl->mcut[inv1];
		T mmax = pimpl->mcut[inv2];

		for (std::size_t i = 0;
			&pimpl->channels[channel].invariants[i] != &invariant; ++i)
		{
			std::size_t const virt = pimpl->channels[channel].invariants[i].in;
			bool const condition = s[virt] > pimpl->mcut[virt] * pimpl->mcut[virt];

			// is there a minimum limit on this invariant?
			if (invariant.lmin.test(i) && condition)
			{
				std::size_t const inv3 = inv1 - virt;
				mmin += sqrt(s[virt]) - pimpl->mcut[inv1] + pimpl->mcut[inv3];
				inv1 = inv3;
			}

			// is there a maximum limit on this invariant?
			if (invariant.lmax.test(i) && condition)
			{
				std::size_t const inv3 = inv2 - virt;
				mmax += sqrt(s[virt]) - pimpl->mcut[inv2] + pimpl->mcut[inv3];
				inv2 = inv3;
			}
		}

		T const smin = mmin * mmin;
		T const smax = (cmf_energy - mmax) * (cmf_energy - mmax);

		// TODO: replace particle_infos calls with model

		s[invariant.in] = map(
			pimpl->particle_infos.at(invariant.idhep).power,
			pimpl->particle_infos.at(invariant.idhep).mass,
			pimpl->particle_infos.at(invariant.idhep).width,
			*r++,
			smin,
			smax
		);
	}

	for (auto const& process : pimpl->channels[channel].processes)
	{
		T const s  = pimpl->s[process.in];
		T const s1 = pimpl->s[process.out1];
		T const s2 = pimpl->s[process.out2];
		T const t1 = pimpl->s[process.in1];
		T const t2 = pimpl->s[process.in2];

		T const lambdas = sqrt(kaellen(s, s1, s2));
		T const lambdat = sqrt(kaellen(s, t1, t2));

		T const tmp = (s + s1 - s2) * (s + t1 - t2);
		T const tmin = s1 + t1 - T(0.5) * (tmp + lambdas * lambdat) / s;
		T       tmax = s1 + t1 - T(0.5) * (tmp - lambdas * lambdat) / s;

		// TODO: make this parameter available from outside
		if (fabs(tmax) < T(1e-7))
		{
			tmax = T();
		}

		T phi = T(2.0) * acos(T(-1.0)) * *r++;
		T const h = map(
			 pimpl->particle_infos.at(process.idhep).power,
			-pimpl->particle_infos.at(process.idhep).mass,
			 pimpl->particle_infos.at(process.idhep).width,
			*r++,
			-tmax,
			-tmin
		);
		T cos_theta = (tmp - T(2.0) * s * (t1 + s1 + h)) / (lambdas * lambdat);

		auto& p1 = p[process.out1];

		T const sqrts = sqrt(s);

		p1 = { T(0.5) * (s + s1 - s2) / sqrts, T(), T(),
			T(0.5) * lambdas / sqrts };

		auto const& q1 = p[process.in1];
		phi = copysign(phi, q1[3]);

		rotate(p1, phi, cos_theta);

		auto const& q2 = p[process.in2];

		std::array<T, 4> const q = {
			q1[0] + q2[0],
			q1[1] + q2[1],
			q1[2] + q2[2],
			q1[3] + q2[3]
		};

		std::array<T, 4> k1 = { q1[0], q1[1], q1[2], q1[3] };

		boost(sqrts, k1, q, false);

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

		T& t = pimpl->s[process.virt];
		t = qt[0] * qt[0] - qt[1] * qt[1] - qt[2] * qt[2] - qt[3] * qt[3];
	}

	for (auto const& decay : pimpl->channels[channel].decays)
	{
		T const s  = pimpl->s[decay.in];
		T const s1 = pimpl->s[decay.out1];
		T const s2 = pimpl->s[decay.out2];

		auto& p1 = p[decay.out1];
		T const sqrts = sqrt(s);

		p1 = { T(0.5) * (s + s1 - s2) / sqrts, T(), T(),
			T(0.5) * sqrt(kaellen(s, s1, s2)) / sqrts };

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

	for (std::size_t i = 2; i != pimpl->particles; ++i)
	{
		momenta.at(4 * i + 0) = p[(1 << i) - 1][0];
		momenta.at(4 * i + 1) = p[(1 << i) - 1][1];
		momenta.at(4 * i + 2) = p[(1 << i) - 1][2];
		momenta.at(4 * i + 3) = p[(1 << i) - 1][3];
	}

	pimpl->cmf_energy_ = cmf_energy;

	// append the extra remaining random numbers to the end of `momenta`
	std::copy_n(r, pimpl->extra_random_numbers,
		std::prev(momenta.end(), pimpl->extra_random_numbers));
}

template <typename T>
luminosity_info<T> lusifer_phase_space_generator<T>::info() const
{
	return luminosity_info<T>(pimpl->cmf_energy_ * pimpl->cmf_energy_);
}

template <typename T>
std::size_t lusifer_phase_space_generator<T>::map_dimensions() const
{
	return 4 * pimpl->particles + pimpl->extra_random_numbers;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class lusifer_constants<double>;
template class lusifer_phase_space_generator<double>;

}
