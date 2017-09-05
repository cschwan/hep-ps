#include "hep/ps/rambo_phase_space_generator.hpp"

#include "catch.hpp"

#include <algorithm>
//#include <iostream>
#include <limits>
#include <random>
#include <vector>

using T = double;

std::array<T, 32> momenta = {
	T( 7.97063578206315810e+01), T( 0.00000000000000000e+00),
	T( 0.00000000000000000e+00), T(-7.97063578206315810e+01),
	T( 7.97063578206315810e+01), T( 0.00000000000000000e+00),
	T( 0.00000000000000000e+00), T( 7.97063578206315810e+01),
	T( 1.40559769420665077e+01), T( 9.87850009899565329e+00),
	T(-6.64796399792236237e+00), T(-7.46929034596835884e+00),
	T( 2.08991451896356182e+01), T(-2.06799352035986743e+01),
	T( 2.05784007638139776e+00), T( 2.20903686082145834e+00),
	T( 9.44895154356087197e+00), T( 1.24401331541535276e+00),
	T( 2.70533374701390494e+00), T( 8.96751278008001407e+00),
	T( 5.48787415399971081e+01), T( 2.54935507955601537e+01),
	T(-3.00343526827748057e+01), T(-3.82059262389317453e+01),
	T( 6.74346435822786550e+00), T(-4.20886431410245798e+00),
	T(-2.36523799876253760e+00), T(-4.70801677411775543e+00),
	T( 5.33864360677752430e+01), T(-1.17272646922700226e+01),
	T( 3.42843808560644092e+01), T( 3.92066837181163876e+01),
};

std::array<T, 1> densities = {
	T(1.07926371985857369e-04)
};

TEST_CASE("numerical momenta", "[rambo_phase_space_generator]")
{
	T const min_energy = T(10.0);
	T const cmf_energy = T(1000.0);

	auto psg = hep::make_rambo_phase_space_generator(
		min_energy,
		cmf_energy,
		6
	);

	REQUIRE( psg->channels()       ==  1 );
	REQUIRE( psg->dimensions()     == 26 );
	REQUIRE( psg->map_dimensions() == 32 );

	std::mt19937 rng;
	std::vector<T> random_numbers(psg->dimensions());

	std::generate(random_numbers.begin(), random_numbers.end(), [&](){
		return std::generate_canonical<T, std::numeric_limits<T>::digits>(rng);
	});

	std::vector<T> p(psg->map_dimensions());
	std::vector<T> d(psg->channels());

//	// generate the test data
//	std::cout.precision(std::numeric_limits<T>::max_digits10);
//	std::cout.setf(std::ios_base::scientific);
//	std::cout.setf(std::ios_base::right);
//
//	std::cout << "std::array<T, " << psg->map_dimensions() << "> momenta = {\n";
//
//	psg->generate(random_numbers, p, 0);
//
//	for (std::size_t i = 0; i != p.size() / 2; ++i)
//	{
//		std::cout << "\tT(" << std::setw(std::cout.precision() + 7)
//			<< p.at(2 * i) << "), T("
//			<< std::setw(std::cout.precision() + 7) << p.at(2 * i + 1)
//			<< "),\n";
//	}
//
//	std::cout << "};\n";
//
//	std::cout << "std::array<T, " << psg->channels() << "> densities = {\n";
//
//	psg->generate(random_numbers, p, 0);
//	psg->densities(d);
//
//	std::cout << "std::array<T, 1> densities = {\n";
//	std::cout << "\tT(" << d.front() << ")\n";
//	std::cout << "};\n";

	psg->generate(random_numbers, p, 0);

	for (std::size_t j = 0; j != momenta.size(); ++j)
	{
		CHECK( p.at(j) == momenta.at(j) );
	}

	psg->densities(d);

	CHECK( d.front() == densities.front() );
}
