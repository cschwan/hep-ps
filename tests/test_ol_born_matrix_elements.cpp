#include "hep/ps/ol_born_matrix_elements.hpp"

#include "catch.hpp"

#include <string>
#include <vector>

using T = HEP_TYPE_T;

TEST_CASE("", "[]")
{
	CHECK_NOTHROW(
		hep::ol_born_matrix_elements<T>(std::vector<std::string>{
			"2 2 -> -11 12 -13 13 1 2",
			"4 4 -> -11 12 -13 13 3 4",
			"2 4 -> -11 12 -13 13 3 2",
			"4 2 -> -11 12 -13 13 1 4"}, 0)
	);
}
