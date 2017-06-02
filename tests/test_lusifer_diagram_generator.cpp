#include "hep/ps/lusifer_constants.hpp"
#include "hep/ps/lusifer_diagram_generator.hpp"

#include "catch.hpp"

#include <string>
#include <vector>

using T = double;

hep::lusifer_constants<T> constants(
	T(125.09), T(4.0e-3),
	T(174.2), T(1.41),
	T(80.385), T(2.085),
	T(91.1876), T(2.4952)
);

TEST_CASE("generator", "[lusifer_diagram_generator]")
{
	auto const diagrams = hep::lusifer_diagram_generator<T>(
		std::vector<std::string>{"sq~uq ne el~nm mu~dq cq~"},
		constants
	);

	CHECK( diagrams.size() == 93 );

	for (auto const& diagram : diagrams)
	{
		CHECK( diagram.propagators().size() == 5 );
		CHECK( diagram.vertices().size() == 6 );
	}
}
