#include "hep/ps/phase_space_generator.hpp"
#include "hep/ps/cofferaa_phase_space_generator.hpp"
#include "hep/ps/lusifer_phase_space_generator.hpp"

#include "catch.hpp"

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

std::size_t histogram(T value, std::vector<std::size_t>& histogram)
{
	using std::fabs;
	using std::pow;

	for (std::size_t power = histogram.size() - 1; power != 0; --power)
	{
		if (fabs(value) < pow(T(10.0), -T(power)))
		{
			++histogram.at(power);

			return power;
		}
	}

	++histogram.front();

	return 0;
}

void run_phase_space_generator(
	hep::phase_space_generator<T>* psg,
	std::vector<std::size_t>& mc_histogram,
	std::vector<std::size_t>& os_histogram,
	std::size_t precision
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

			CAPTURE( i );
			CAPTURE( channel );

			// pick the largest absolute value
			T const value = *std::max_element(sums.begin(), sums.end(),
				[](T one, T two) { return fabs(one) < fabs(two); });

			// check momentum conservation
			CHECK( histogram(value, mc_histogram) > precision );

			for (std::size_t particle = 0; particle != particles; ++particle)
			{
				T const invariant =
					p.at(4 * particle + 0) * p.at(4 * particle + 0) -
					p.at(4 * particle + 1) * p.at(4 * particle + 1) -
					p.at(4 * particle + 2) * p.at(4 * particle + 2) -
					p.at(4 * particle + 3) * p.at(4 * particle + 3);

				CAPTURE( particle );

				// check on-shellness
				CHECK( histogram(invariant, os_histogram) > precision );
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

	run_phase_space_generator(psg.get(), mc_histogram, os_histogram, 9);

	CHECK( mc_histogram.at( 0) == 0 );
	CHECK( mc_histogram.at( 1) == 0 );
	CHECK( mc_histogram.at( 2) == 0 );
	CHECK( mc_histogram.at( 3) == 0 );
	CHECK( mc_histogram.at( 4) == 0 );
	CHECK( mc_histogram.at( 5) == 0 );
	CHECK( mc_histogram.at( 6) == 0 );
	CHECK( mc_histogram.at( 7) == 0 );
	CHECK( mc_histogram.at( 8) == 0 );
	CHECK( mc_histogram.at( 9) == 0 );
	CHECK( mc_histogram.at(10) == 0 );
	CHECK( mc_histogram.at(11) == 0 );
	CHECK( mc_histogram.at(12) == 18087 );
	CHECK( mc_histogram.at(13) == 184728 );
	CHECK( mc_histogram.at(14) == 206748 );
	CHECK( mc_histogram.at(15) == 45555 );
	CHECK( mc_histogram.at(16) == 127 );
	CHECK( mc_histogram.at(17) == 0 );
	CHECK( mc_histogram.at(18) == 0 );
	CHECK( mc_histogram.at(19) == 0 );
	CHECK( mc_histogram.at(20) == 0 );
	CHECK( mc_histogram.at(21) == 0 );
	CHECK( mc_histogram.at(22) == 0 );
	CHECK( mc_histogram.at(23) == 0 );
	CHECK( mc_histogram.at(24) == 0 );
	CHECK( mc_histogram.at(25) == 0 );
	CHECK( mc_histogram.at(26) == 0 );
	CHECK( mc_histogram.at(27) == 0 );
	CHECK( mc_histogram.at(28) == 0 );
	CHECK( mc_histogram.at(29) == 755 );

	CHECK( os_histogram.at( 0) == 0 );
	CHECK( os_histogram.at( 1) == 0 );
	CHECK( os_histogram.at( 2) == 0 );
	CHECK( os_histogram.at( 3) == 0 );
	CHECK( os_histogram.at( 4) == 0 );
	CHECK( os_histogram.at( 5) == 0 );
	CHECK( os_histogram.at( 6) == 0 );
	CHECK( os_histogram.at( 7) == 0 );
	CHECK( os_histogram.at( 8) == 0 );
	CHECK( os_histogram.at( 9) == 0 );
	CHECK( os_histogram.at(10) == 25071 );
	CHECK( os_histogram.at(11) == 180500 );
	CHECK( os_histogram.at(12) == 573990 );
	CHECK( os_histogram.at(13) == 432003 );
	CHECK( os_histogram.at(14) == 406643 );
	CHECK( os_histogram.at(15) == 368135 );
	CHECK( os_histogram.at(16) == 240630 );
	CHECK( os_histogram.at(17) == 152021 );
	CHECK( os_histogram.at(18) == 114142 );
	CHECK( os_histogram.at(19) == 62151 );
	CHECK( os_histogram.at(20) == 33942 );
	CHECK( os_histogram.at(21) == 20295 );
	CHECK( os_histogram.at(22) == 8343 );
	CHECK( os_histogram.at(23) == 4387 );
	CHECK( os_histogram.at(24) == 2081 );
	CHECK( os_histogram.at(25) == 864 );
	CHECK( os_histogram.at(26) == 177 );
	CHECK( os_histogram.at(27) == 13 );
	CHECK( os_histogram.at(28) == 0 );
	CHECK( os_histogram.at(29) == 1478612 );
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

	run_phase_space_generator(psg.get(), mc_histogram, os_histogram, 8);

	CHECK( mc_histogram.at( 0) == 0 );
	CHECK( mc_histogram.at( 1) == 0 );
	CHECK( mc_histogram.at( 2) == 0 );
	CHECK( mc_histogram.at( 3) == 0 );
	CHECK( mc_histogram.at( 4) == 0 );
	CHECK( mc_histogram.at( 5) == 0 );
	CHECK( mc_histogram.at( 6) == 0 );
	CHECK( mc_histogram.at( 7) == 0 );
	CHECK( mc_histogram.at( 8) == 0 );
	CHECK( mc_histogram.at( 9) == 0 );
	CHECK( mc_histogram.at(10) == 0 );
	CHECK( mc_histogram.at(11) == 0 );
	CHECK( mc_histogram.at(12) == 1064 );
	CHECK( mc_histogram.at(13) == 31203 );
	CHECK( mc_histogram.at(14) == 43341 );
	CHECK( mc_histogram.at(15) == 16633 );
	CHECK( mc_histogram.at(16) == 124 );
	CHECK( mc_histogram.at(17) == 4 );
	CHECK( mc_histogram.at(18) == 0 );
	CHECK( mc_histogram.at(19) == 0 );
	CHECK( mc_histogram.at(20) == 0 );
	CHECK( mc_histogram.at(21) == 0 );
	CHECK( mc_histogram.at(22) == 0 );
	CHECK( mc_histogram.at(23) == 0 );
	CHECK( mc_histogram.at(24) == 0 );
	CHECK( mc_histogram.at(25) == 0 );
	CHECK( mc_histogram.at(26) == 0 );
	CHECK( mc_histogram.at(27) == 0 );
	CHECK( mc_histogram.at(28) == 0 );
	CHECK( mc_histogram.at(29) == 631 );

	CHECK( os_histogram.at( 0) == 0 );
	CHECK( os_histogram.at( 1) == 0 );
	CHECK( os_histogram.at( 2) == 0 );
	CHECK( os_histogram.at( 3) == 0 );
	CHECK( os_histogram.at( 4) == 0 );
	CHECK( os_histogram.at( 5) == 0 );
	CHECK( os_histogram.at( 6) == 0 );
	CHECK( os_histogram.at( 7) == 0 );
	CHECK( os_histogram.at( 8) == 0 );
	CHECK( os_histogram.at( 9) == 302 );
	CHECK( os_histogram.at(10) == 21613 );
	CHECK( os_histogram.at(11) == 85790 );
	CHECK( os_histogram.at(12) == 117976 );
	CHECK( os_histogram.at(13) == 96209 );
	CHECK( os_histogram.at(14) == 85897 );
	CHECK( os_histogram.at(15) == 53988 );
	CHECK( os_histogram.at(16) == 28717 );
	CHECK( os_histogram.at(17) == 15063 );
	CHECK( os_histogram.at(18) == 9558 );
	CHECK( os_histogram.at(19) == 4667 );
	CHECK( os_histogram.at(20) == 2088 );
	CHECK( os_histogram.at(21) == 977 );
	CHECK( os_histogram.at(22) == 411 );
	CHECK( os_histogram.at(23) == 191 );
	CHECK( os_histogram.at(24) == 70 );
	CHECK( os_histogram.at(25) == 18 );
	CHECK( os_histogram.at(26) == 0 );
	CHECK( os_histogram.at(27) == 0 );
	CHECK( os_histogram.at(28) == 0 );
	CHECK( os_histogram.at(29) == 220465 );
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

	run_phase_space_generator(psg.get(), mc_histogram, os_histogram, 4);

	CHECK( mc_histogram.at( 0) == 0 );
	CHECK( mc_histogram.at( 1) == 0 );
	CHECK( mc_histogram.at( 2) == 0 );
	CHECK( mc_histogram.at( 3) == 0 );
	CHECK( mc_histogram.at( 4) == 0 );
	CHECK( mc_histogram.at( 5) == 0 );
	CHECK( mc_histogram.at( 6) == 0 );
	CHECK( mc_histogram.at( 7) == 0 );
	CHECK( mc_histogram.at( 8) == 0 );
	CHECK( mc_histogram.at( 9) == 0 );
	CHECK( mc_histogram.at(10) == 592 );
	CHECK( mc_histogram.at(11) == 15143 );
	CHECK( mc_histogram.at(12) == 23380 );
	CHECK( mc_histogram.at(13) == 21547 );
	CHECK( mc_histogram.at(14) == 22780 );
	CHECK( mc_histogram.at(15) == 8020 );
	CHECK( mc_histogram.at(16) == 86 );
	CHECK( mc_histogram.at(17) == 0 );
	CHECK( mc_histogram.at(18) == 0 );
	CHECK( mc_histogram.at(19) == 0 );
	CHECK( mc_histogram.at(20) == 0 );
	CHECK( mc_histogram.at(21) == 0 );
	CHECK( mc_histogram.at(22) == 0 );
	CHECK( mc_histogram.at(23) == 0 );
	CHECK( mc_histogram.at(24) == 0 );
	CHECK( mc_histogram.at(25) == 0 );
	CHECK( mc_histogram.at(26) == 0 );
	CHECK( mc_histogram.at(27) == 0 );
	CHECK( mc_histogram.at(28) == 0 );
	CHECK( mc_histogram.at(29) == 1452 );

	CHECK( os_histogram.at( 0) == 0 );
	CHECK( os_histogram.at( 1) == 0 );
	CHECK( os_histogram.at( 2) == 0 );
	CHECK( os_histogram.at( 3) == 0 );
	CHECK( os_histogram.at( 4) == 0 );
	CHECK( os_histogram.at( 5) == 246 );
	CHECK( os_histogram.at( 6) == 12335 );
	CHECK( os_histogram.at( 7) == 34942 );
	CHECK( os_histogram.at( 8) == 46704 );
	CHECK( os_histogram.at( 9) == 59669 );
	CHECK( os_histogram.at(10) == 59129 );
	CHECK( os_histogram.at(11) == 69356 );
	CHECK( os_histogram.at(12) == 68524 );
	CHECK( os_histogram.at(13) == 51682 );
	CHECK( os_histogram.at(14) == 44265 );
	CHECK( os_histogram.at(15) == 27381 );
	CHECK( os_histogram.at(16) == 14747 );
	CHECK( os_histogram.at(17) == 8218 );
	CHECK( os_histogram.at(18) == 5365 );
	CHECK( os_histogram.at(19) == 2416 );
	CHECK( os_histogram.at(20) == 1200 );
	CHECK( os_histogram.at(21) == 426 );
	CHECK( os_histogram.at(22) == 286 );
	CHECK( os_histogram.at(23) == 104 );
	CHECK( os_histogram.at(24) == 23 );
	CHECK( os_histogram.at(25) == 0 );
	CHECK( os_histogram.at(26) == 0 );
	CHECK( os_histogram.at(27) == 0 );
	CHECK( os_histogram.at(28) == 0 );
	CHECK( os_histogram.at(29) == 236982 );
}
