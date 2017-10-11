#include "hep/ps/parton.hpp"
#include "hep/ps/parton_dfs.hpp"

#include "test_structures.hpp"

#include "catch.hpp"

#include <random>
#include <vector>

using T = HEP_TYPE_T;

TEST_CASE("parton distribution functions", "[parton_dfs]")
{
	hep::parton_dfs<T> pdfs{"CT14nlo"};

	T const scale = T(100.0);

	std::mt19937 rng;
	std::uniform_real_distribution<T> dist{T(), T(1.0)};
	std::vector<hep::parton_array<T>> results1(pdfs.count());

	for (std::size_t i = 0; i != 100; ++i)
	{
		T const x = dist(rng);

		// this evaluation call should store the value of all PDFs in
		// `results1`, the first element giving the central PDF
		pdfs.eval(x, scale, results1);

		// this evaluation call should return only the central PDF
		auto const result2 = pdfs.eval(x, scale);

		// check if our assumptions are true
		for (auto const parton : hep::parton_list())
		{
			CHECK( results1.front()[parton] == result2[parton] );
		}
	}
}
