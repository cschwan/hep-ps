#include "hep/ps/fortran_helper.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "lusifer_interfaces.hpp"

#include <cassert>
#include <cmath>

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
	std::size_t channels;
	std::size_t particles;
	std::size_t max_particles;
	T cmf_energy_;
};

template <typename T>
lusifer_phase_space_generator<T>::lusifer_phase_space_generator(
	std::string const& process,
	lusifer_constants<T> const& constants
)
	: pimpl(new impl())
{
	// each particle must be specified with three characters
	assert( process.size() % 3 == 0 );

	pimpl->particles = process.size() / 3;

	// there must be at least four particles
	assert( pimpl->particles >= 4 );

	int maxex;
	int maxgen;
	lusifer_extra_generatormax(&maxex, &maxgen);

	pimpl->max_particles = maxex;

	// the number of particles must be supported by the generator
	assert( pimpl->max_particles >= pimpl->particles );

	// FORTRAN counting: use the first generator
	int generator = 1;
	int nex = pimpl->particles;
	double mw = constants.mass_w;
	double gw = constants.width_w;
	double mz = constants.mass_z;
	double gz = constants.width_z;
	double mh = constants.mass_h;
	double gh = constants.width_h;
	double mt = constants.mass_t;
	double gt = constants.width_t;

	// set constants and the number of particles
	lusifer_extra_generatorset(
		&generator,
		&nex,
		&mw,
		&gw,
		&mz,
		&gz,
		&mh,
		&gh,
		&mt,
		&gt
	);

	// fill up the string with three spaces for particle not used
	std::string process0 = process;
	process0.append(' ', 3 * (maxex - pimpl->particles));

	// TODO: what is the meaning of this parameter?
	int lightfermions = 0;
	// do not include cuts in the phase space generation
	int includecuts = 0;
	// do not print channel information
	int sout = 0;

	lusifer_initphasespace(
		process0.c_str(),
		&generator,
		&lightfermions,
		&includecuts,
		&sout,
		process0.size()
	);

	int channels;
	lusifer_extra_generatordata(&generator, &channels);

	// there must be at least one channel, otherwise something went wrong
	assert( channels > 0 );

	pimpl->channels = channels;
}

template <typename T>
lusifer_phase_space_generator<T>::~lusifer_phase_space_generator() = default;

template <typename T>
std::size_t lusifer_phase_space_generator<T>::channels() const
{
	return pimpl->channels;
}

template <typename T>
T lusifer_phase_space_generator<T>::densities(std::vector<T>& densities)
{
	assert( densities.size() >= channels() );

	int generator = 1;
	int switch_ = 2;

	lusifer_density(densities.data(), &generator, &switch_);

	using std::acos;
	using std::pow;

	return pow(T(0.5) / acos(T(-1.0)), T(3 * pimpl->particles - 10));
}

template <typename T>
std::size_t lusifer_phase_space_generator<T>::dimensions() const
{
	return 3 * (pimpl->particles - 4) + 2;
}

template <typename T>
void lusifer_phase_space_generator<T>::generate(
	std::vector<T> const& random_numbers,
	std::vector<T>& momenta,
	T cmf_energy,
	std::size_t channel
) {
	assert( random_numbers.size() >= dimensions() );
	assert( momenta.size() == map_dimensions() );

	double const value = 0.5 * cmf_energy;
	// incoming beam momenta
	double kbeam[] = { value, value, 0.0, 0.0, 0.0, 0.0, -value, value };
	double x1 = 1.0;
	double x2 = 1.0;
	// dummy parameter
	double g = 0.0;
	// FORTRAN counting again
	int channel_ = channel + 1;
	int generator = 1;
	int switch_ = 1;

	// the vector containing the random numbers must have the maximimum size
	std::vector<T> random_numbers0;
	random_numbers0.reserve(3 * (pimpl->max_particles - 4) + 2);
	random_numbers0 = random_numbers;
	random_numbers0.resize(3 * (pimpl->max_particles - 4) + 2);

	// the same holds true for the momenta
	std::vector<T> momenta0;
	momenta0.reserve(4 * pimpl->max_particles);
	momenta0 = momenta;
	momenta0.resize(4 * pimpl->max_particles);

	lusifer_phasespace(
		random_numbers0.data(),
		kbeam,
		momenta0.data(),
		&x1,
		&x2,
		&g,
		&channel_,
		&generator,
		&switch_
	);

	fortran_ordering_to_cpp(momenta0);
	momenta0.resize(momenta.size());
	momenta = momenta0;

	pimpl->cmf_energy_ = cmf_energy;
}

template <typename T>
luminosity_info<T> lusifer_phase_space_generator<T>::info() const
{
	return luminosity_info<T>(pimpl->cmf_energy_ * pimpl->cmf_energy_);
}

template <typename T>
std::size_t lusifer_phase_space_generator<T>::map_dimensions() const
{
	return 4 * pimpl->particles;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class lusifer_constants<double>;
template class lusifer_phase_space_generator<double>;

}
