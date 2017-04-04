#include "hep/ps/cofferaa_phase_space_generator.hpp"
#include "hep/ps/fortran_helper.hpp"

#include "cofferaa_interfaces.hpp"

#include <cassert>

namespace hep
{

template <typename T>
class cofferaa_phase_space_generator<T>::impl
{
public:
	T energy;

	std::size_t channels;
	std::size_t max_particles;
	std::size_t particles;

	std::vector<T> current_point;
};

template <typename T>
cofferaa_phase_space_generator<T>::cofferaa_phase_space_generator(
	std::vector<int> const& process,
	lusifer_constants<T> const& constants,
	std::vector<std::tuple<int, int, int>> const& dipoles,
	std::size_t extra_random_numbers
)
	: pimpl(new impl())
{
	// TODO: NYI
	assert( extra_random_numbers == 0 );

	pimpl->particles = process.size();

	int maxex;
	int maxgen;
	cofferaa_extra_max(&maxex, &maxgen);

	pimpl->max_particles = maxex;

	assert( pimpl->max_particles >= pimpl->particles );

	pimpl->current_point.resize(4 * maxex);

	// FORTRAN counting: use the first generator
	int generator = 1;
	double mw = constants.mass_w;
	double gw = constants.width_w;
	double mz = constants.mass_z;
	double gz = constants.width_z;
	double mh = constants.mass_h;
	double gh = constants.width_h;
	double mt = constants.mass_t;
	double gt = constants.width_t;

	// set constants and the number of particles
	cofferaa_extra_set(
		&generator,
		&mw,
		&gw,
		&mz,
		&gz,
		&mh,
		&gh,
		&mt,
		&gt
	);

	// only used for technical cuts
	double energy = 0.0;
	// output variable; minimum energy that the generator produces
	double smin;
	// number of external particles
	int next = process.size();
	// EW SM + QCD
	int smodel = 12;
	int sincludecuts = 0;
	// generate phase space according to the dipole mappings
	int ssub = (dipoles.size() != 0) ? 1 : 0;

	std::vector<int> dipole_emitter;
	std::vector<int> dipole_unresolved;
	std::vector<int> dipole_spectator;

	for (auto const& dipole : dipoles)
	{
		dipole_emitter.push_back(std::get<0>(dipole));
		dipole_unresolved.push_back(std::get<1>(dipole));
		dipole_spectator.push_back(std::get<2>(dipole));
	}

	int dipole_count = dipoles.size();

	cofferaa_initgenerator(
		&energy,
		&smin,
		process.data(),
		&generator,
		&next,
		&smodel,
		&sincludecuts,
		&ssub,
		&dipole_count,
		dipole_emitter.data(),
		dipole_unresolved.data(),
		dipole_spectator.data()
	);

	int channels;
	cofferaa_extra_data(&generator, &channels);

	assert( channels > 0 );

	pimpl->channels = channels;
}

template <typename T>
cofferaa_phase_space_generator<T>::cofferaa_phase_space_generator(
	cofferaa_phase_space_generator&& psg
)
	: pimpl(std::move(psg.pimpl))
{
}

template <typename T>
cofferaa_phase_space_generator<T>::~cofferaa_phase_space_generator() = default;

template <typename T>
std::size_t cofferaa_phase_space_generator<T>::channels() const
{
	return pimpl->channels;
}

template <typename T>
T cofferaa_phase_space_generator<T>::densities(std::vector<T>& densities)
{
	assert( densities.size() >= channels() );

	int generator = 1;
	int switch_ = 2;

	cofferaa_density(
		pimpl->current_point.data(),
		densities.data(),
		&generator,
		&switch_
	);

	using std::acos;
	using std::pow;

	return pow(T(0.5) / acos(T(-1.0)), T(3 * pimpl->particles - 10));
}

template <typename T>
std::size_t cofferaa_phase_space_generator<T>::dimensions() const
{
	return 3 * (pimpl->particles - 4) + 2;
}

template <typename T>
void cofferaa_phase_space_generator<T>::generate(
	std::vector<T> const& random_numbers,
	std::vector<T>& momenta,
	T cmf_energy,
	std::size_t channel
) {
	assert( random_numbers.size() >= dimensions() );
	assert( momenta.size() == map_dimensions() );

	double const value = 0.5 * static_cast <double> (cmf_energy);
	double kbeam[] = { value, value, 0.0, 0.0, 0.0, 0.0, -value, value };
	// dummy parameter
	double g;
	// FORTRAN counting
	int channel_ = channel + 1;
	int generator = 1;
	int switch_ = 1;

	cofferaa_generation(
		random_numbers.data(),
		&kbeam[0],
		pimpl->current_point.data(),
		&g,
		&channel_,
		&generator,
		&switch_
	);

	momenta = pimpl->current_point;
	fortran_ordering_to_cpp(momenta);
	momenta.resize(map_dimensions());
}

template <typename T>
luminosity_info<T> cofferaa_phase_space_generator<T>::info() const
{
	return luminosity_info<T>(pimpl->energy * pimpl->energy);
}

template <typename T>
std::size_t cofferaa_phase_space_generator<T>::map_dimensions() const
{
	return 4 * pimpl->particles;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class cofferaa_phase_space_generator<double>;

}