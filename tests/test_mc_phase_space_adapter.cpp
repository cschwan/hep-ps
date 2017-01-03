#include "hep/mc/multi_channel.hpp"
#include "hep/mc/multi_channel_integrand.hpp"
#include "hep/mc/multi_channel_point.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "hep/ps/mc_phase_space_adapter.hpp"

#include "catch.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

using T = double;

hep::lusifer_constants<T> constants(
	T(125.09), T(4.0e-3),
	T(174.2), T(1.41),
	T(80.385), T(2.085),
	T(91.1876), T(2.4952)
);

template <typename T>
T s_invariant(std::vector<T> const& momenta, std::size_t i, std::size_t j)
{
	std::array<T, 4> const vector = {
		momenta.at(4 * i + 0) + momenta.at(4 * j + 0),
		momenta.at(4 * i + 1) + momenta.at(4 * j + 1),
		momenta.at(4 * i + 2) + momenta.at(4 * j + 2),
		momenta.at(4 * i + 3) + momenta.at(4 * j + 3)
	};

	T const invariant =
		vector.at(0) * vector.at(0) -
		vector.at(1) * vector.at(1) -
		vector.at(2) * vector.at(2) -
		vector.at(3) * vector.at(3);

	return invariant;
}

T t_invariant(std::vector<T> const& momenta, std::size_t i, std::size_t j)
{
	std::array<T, 4> const vector = {
		momenta.at(4 * i + 0) - momenta.at(4 * j + 0),
		momenta.at(4 * i + 1) - momenta.at(4 * j + 1),
		momenta.at(4 * i + 2) - momenta.at(4 * j + 2),
		momenta.at(4 * i + 3) - momenta.at(4 * j + 3)
	};

	T const invariant =
		vector.at(0) * vector.at(0) -
		vector.at(1) * vector.at(1) -
		vector.at(2) * vector.at(2) -
		vector.at(3) * vector.at(3);

	return invariant;
}

TEST_CASE("muon pair creation", "[mc_phase_space_adapter]")
{
	using std::acos;

	T const pi = acos(T(-1.0));
	T const alpha = T(1.0 / 137.035999074);
	T const hbarc2 = T(0.3893793656e+6);

	auto function = [=](hep::multi_channel_point<T> const& x) {
		CHECK( x.channel() < 2 );
		CHECK( x.coordinates().size() == 16 );
		CHECK( x.point().size() == 2 );
		CHECK( x.weight() >= T() );

		T const s = s_invariant(x.coordinates(), 0, 1);
		T const t = t_invariant(x.coordinates(), 0, 2);
		T const u = t_invariant(x.coordinates(), 0, 3);

		// matrix element squared with correct averaging factor
		T const me2 = T(32.0) * pi*pi * alpha*alpha * (t*t + u*u) / (s*s);
		// flux factor
		T const flux_factor = T(0.5) / s;

		return hbarc2 * flux_factor * me2;
	};

	// 1000 GeV = 1 TeV
	T const cmf_energy = T(1000.0);

	auto generator = hep::mc_phase_space_adapter<
		hep::lusifer_phase_space_generator<T>>(
			cmf_energy,
			"el~el mu mu~",
			constants
	);

	auto const results = hep::multi_channel(
		hep::make_multi_channel_integrand<T>(
			function,
			generator.dimensions(),
			std::ref(generator),
			generator.map_dimensions(),
			generator.channels()
		),
		std::vector<std::size_t>(2, 1000)
	);

	// result is 1 R, with R = 4*pi*alpha^2/(3*s) where we set alpha=1
	T const oneR = T(4.0) * pi * alpha * alpha * hbarc2 /
		(T(3.0) * cmf_energy * cmf_energy);

	T const value = results.back().value();
	T const error = results.back().error();

	CHECK( (value / oneR) == Approx(T(1.0087047847457573)) );
}
