#include "hep/ps/phase_space_generator.hpp"
#include "hep/ps/cofferaa_phase_space_generator.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"

#include <catch.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

using T = double;

hep::lusifer_constants<T> constants(
	T(125.09), T(4.0e-3),
	T(174.2), T(1.41),
	T(80.385), T(2.085),
	T(91.1876), T(2.4952)
);

void histogram(T value, std::vector<std::size_t>& histogram)
{
	using std::fabs;
	using std::pow;

	for (std::size_t power = histogram.size() - 1; power != 0; --power)
	{
		if (fabs(value) < pow(T(10.0), -T(power)))
		{
			++histogram.at(power);

			return;
		}
	}

	++histogram.front();
}

void run_phase_space_generator(
	hep::phase_space_generator<T>* psg,
	std::vector<std::size_t>& mc_histogram,
	std::vector<std::size_t>& os_histogram
) {
	std::mt19937 rng;
	std::vector<T> random_numbers(psg->dimensions());

	std::vector<T> p(psg->map_dimensions());
	std::size_t const particles = psg->map_dimensions() / 4;

	for (std::size_t i = 0; i != 1000; ++i)
	{
		std::generate(random_numbers.begin(), random_numbers.end(), [&](){
			return std::generate_canonical<T,
				std::numeric_limits<T>::digits>(rng);
		});

		for (std::size_t channel = 0; channel != psg->channels(); ++channel)
		{
			psg->generate(random_numbers, p, channel);

			std::array<T, 4> sums = { T(), T(), T(), T() };

			for (std::size_t particle = 0; particle != particles; ++particle)
			{
				// count incoming particles as negative
				T const sign = (particle < 2) ? T(-1.0) : T(1.0);

				sums.at(0) += sign * p.at(4 * particle + 0);
				sums.at(1) += sign * p.at(4 * particle + 1);
				sums.at(2) += sign * p.at(4 * particle + 2);
				sums.at(3) += sign * p.at(4 * particle + 3);
			}

			// pick the largest absolute value
			T const value = *std::max_element(sums.begin(), sums.end(),
				[](T one, T two) { return fabs(one) < fabs(two); });

			// check momentum conservation
			histogram(value, mc_histogram);

			for (std::size_t particle = 0; particle != particles; ++particle)
			{
				T const invariant =
					p.at(4 * particle + 0) * p.at(4 * particle + 0) -
					p.at(4 * particle + 1) * p.at(4 * particle + 1) -
					p.at(4 * particle + 2) * p.at(4 * particle + 2) -
					p.at(4 * particle + 3) * p.at(4 * particle + 3);

				histogram(invariant, os_histogram);
			}
		}
	}
}

TEST_CASE("Cofferaa 1 TeV VBS Real QCD PS Quality", "[]")
{
	T const min_energy = T(10.0);
	T const cmf_energy = T(1000.0);

	auto psg = hep::make_cofferaa_phase_space_generator(
		min_energy,
		cmf_energy,
		std::vector<int>{-3, 2, 12, -11, 14, -13, 1, -4, 26},
		constants
	);

	std::vector<std::size_t> mc_histogram(30);
	std::vector<std::size_t> os_histogram(30);

	std::vector<std::size_t> reference_mc = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18087, 184728, 206748, 45555, 127,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 755
	};
	std::vector<std::size_t> reference_os = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25071, 180500, 573990, 432003, 406643,
		368135, 240630, 152021, 114142, 62151, 33942, 20295, 8343, 4387, 2081,
		864, 177, 13, 0, 1478612
	};

	run_phase_space_generator(psg.get(), mc_histogram, os_histogram);

	CHECK_THAT( mc_histogram , Catch::Equals(reference_mc) );
	CHECK_THAT( os_histogram , Catch::Equals(reference_os) );
}

TEST_CASE("Lusifer 1 TeV VBS Real QCD PS Quality", "[]")
{
	T const min_energy = T(10.0);
	T const cmf_energy = T(1000.0);

	auto psg = hep::make_lusifer_phase_space_generator(
		min_energy,
		cmf_energy,
		"sq~uq ne el~nm mu~dq cq~",
		constants
	);

	std::vector<std::size_t> mc_histogram(30);
	std::vector<std::size_t> os_histogram(30);

	std::vector<std::size_t> reference_mc = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1066, 31366, 43093, 16695, 148, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 632
	};
	std::vector<std::size_t> reference_os = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 328, 21062, 85918, 118256, 95875, 86583,
		54428, 28192, 15189, 9330, 4655, 2021, 1014, 434, 190, 79, 22, 0, 0, 0,
		220424
	};

	run_phase_space_generator(psg.get(), mc_histogram, os_histogram);

	CHECK_THAT( mc_histogram , Catch::Equals(reference_mc) );
	CHECK_THAT( os_histogram , Catch::Equals(reference_os) );
}

TEST_CASE("Lusifer 100 TeV VBS Real QCD PS Quality", "[]")
{
	T const min_energy = T(10.0);
	T const cmf_energy = T(100000.0);

	auto psg = hep::make_lusifer_phase_space_generator(
		min_energy,
		cmf_energy,
		"sq~uq ne el~nm mu~dq cq~",
		constants
	);

	std::vector<std::size_t> mc_histogram(30);
	std::vector<std::size_t> os_histogram(30);

	std::vector<std::size_t> reference_mc = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 577, 15142, 23494, 21386, 22932, 8031,
		65, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1373
	};
	std::vector<std::size_t> reference_os = {
		0, 0, 0, 0, 0, 249, 12600, 34747, 46823, 59468, 58567, 69558, 68817,
		51616, 43979, 27374, 15095, 8059, 5160, 2502, 1224, 426, 257, 105, 15,
		0, 0, 0, 0, 237359
	};

	run_phase_space_generator(psg.get(), mc_histogram, os_histogram);

	CHECK_THAT( mc_histogram , Catch::Equals(reference_mc) );
	CHECK_THAT( os_histogram , Catch::Equals(reference_os) );
}
