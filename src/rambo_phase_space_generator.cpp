#include "hep/ps/boost.hpp"
#include "hep/ps/rambo_phase_space_generator.hpp"

#include "hadron_hadron_psg_adapter.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iterator>
#include <vector>

namespace
{

template <typename T>
class rambo_psg
{
public:
	using numeric_type = T;

	rambo_psg(
		std::size_t final_state_particles,
		std::size_t extra_random_numbers
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
	std::size_t final_state_particles_;
	std::size_t extra_random_numbers_;
	T cmf_energy_;
};

template <typename T>
rambo_psg<T>::rambo_psg(
	std::size_t final_state_particles,
	std::size_t extra_random_numbers
)
	: final_state_particles_{final_state_particles}
	, extra_random_numbers_{extra_random_numbers}
{
}

template <typename T>
std::size_t rambo_psg<T>::channels() const
{
	return 1;
}

template <typename T>
T rambo_psg<T>::densities(std::vector<T>& densities)
{
	using std::acos;
	using std::pow;
	using std::tgamma;

	std::size_t const n = final_state_particles_;
	T const w = cmf_energy_;

	densities.front()  = T(8.0) * acos(T(-1.0));
	densities.front() *= (T(n - 1) * tgamma(T(n - 1)) * tgamma(T(n - 1)));
	densities.front() /= pow(T(0.25) * w / acos(T(-1.0)), T(2 * n - 4));

	return T(1.0);
}

template <typename T>
std::size_t rambo_psg<T>::dimensions() const
{
	return 4 * final_state_particles_ + extra_random_numbers_;
}

template <typename T>
void rambo_psg<T>::generate(
	std::vector<T> const& random_numbers,
	std::vector<T>& momenta,
	T cmf_energy,
	std::size_t channel
) {
	assert( channel == 0 );

	using std::acos;
	using std::cos;
	using std::log;
	using std::sin;
	using std::sqrt;

	cmf_energy_ = cmf_energy;

	T const w = cmf_energy;
	std::size_t const n = final_state_particles_;

	momenta.at(0) = T(0.5) * w;
	momenta.at(1) = T{};
	momenta.at(2) = T{};
	momenta.at(3) = T(-0.5) * w;
	momenta.at(4) = T(0.5) * w;
	momenta.at(5) = T{};
	momenta.at(6) = T{};
	momenta.at(7) = T(0.5) * w;

	std::array<T, 4> Q{ T{}, T{}, T{}, T{} };

	// generate isotropic momenta (equation 3.1)
	for (std::size_t i = 0; i != n; ++i)
	{
		// extract random numbers
		T const rho0 = random_numbers.at(4 * i + 0);
		T const rho1 = random_numbers.at(4 * i + 1);
		T const rho2 = random_numbers.at(4 * i + 2);
		T const rho3 = random_numbers.at(4 * i + 3);

		// get a random cosine of an azimuthal angle between 0 and pi
		T const cos_theta = T(2.0) * rho0 - T(1.0);
		// ... compute the corresponding sine (always positive)
		T const sin_theta = sqrt(T(1.0) - cos_theta * cos_theta);
		// get a random polar angle
		T const phi = T(2.0) * acos(T(-1.0)) * rho1;
		// get a random length for the vector
		T const length = -log(rho2 * rho3);

		// assemble the components
		momenta.at(4 * (i + 2) + 0) = length;
		momenta.at(4 * (i + 2) + 1) = length * sin_theta * cos(phi);
		momenta.at(4 * (i + 2) + 2) = length * sin_theta * sin(phi);
		momenta.at(4 * (i + 2) + 3) = length * cos_theta;

		Q[0] += momenta.at(4 * (i + 2) + 0);
		Q[1] += momenta.at(4 * (i + 2) + 1);
		Q[2] += momenta.at(4 * (i + 2) + 2);
		Q[3] += momenta.at(4 * (i + 2) + 3);
	}

	T const M = sqrt(Q[0] * Q[0] - Q[1] * Q[1] - Q[2] * Q[2] - Q[3] * Q[3]);
	T const x = w / M;

	for (std::size_t i = 0; i != n; ++i)
	{
		std::array<T, 4> p = {
			momenta.at(4 * (i + 2) + 0),
			momenta.at(4 * (i + 2) + 1),
			momenta.at(4 * (i + 2) + 2),
			momenta.at(4 * (i + 2) + 3)
		};

		// boost the momenta back into the center-of-mass frame
		hep::boost(M, p, Q, false);

		momenta.at(4 * (i + 2) + 0) = x * p[0];
		momenta.at(4 * (i + 2) + 1) = x * p[1];
		momenta.at(4 * (i + 2) + 2) = x * p[2];
		momenta.at(4 * (i + 2) + 3) = x * p[3];
	}

	std::copy_n(
		std::prev(random_numbers.end(), extra_random_numbers_),
		extra_random_numbers_,
		std::prev(momenta.end(), extra_random_numbers_)
	);
}

template <typename T>
std::size_t rambo_psg<T>::map_dimensions() const
{
	return 4 * (final_state_particles_ + 2) + extra_random_numbers_;
}

}

namespace hep
{

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_rambo_phase_space_generator(
	T min_energy,
	T cmf_energy,
	std::size_t final_state_particles,
	std::size_t extra_random_numbers
) {
	return std::make_unique<hadron_hadron_psg_adapter<rambo_psg<T>>>(
		min_energy,
		cmf_energy,
		final_state_particles,
		extra_random_numbers
	);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template std::unique_ptr<phase_space_generator<double>>
make_rambo_phase_space_generator(
	double,
	double,
	std::size_t,
	std::size_t
);

}
