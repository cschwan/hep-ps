#include "hep/ps/tensor_decompose.hpp"

#include <array>
#include <cstddef>

#include "catch.hpp"

using T = HEP_TYPE_T;

T lhs_matrix_component(
	std::array<T, 4> const& p,
	T denominator,
	std::size_t i,
	std::size_t j
) {
	T const factor = denominator * (p[0] * p[0] - p[1] * p[1] - p[2] * p[2] -
		p[3] * p[3]);

	T result = p[i] * p[j] / factor;

	if (i == j)
	{
		result += (i == 0) ? T(-1.0) : T(1.0);
	}

	return result;
}

T rhs_matrix_component(
	std::array<std::array<T, 4>, 4> const& results,
	std::array<T, 4> const& factors,
	std::size_t i,
	std::size_t j
) {
	T result = T();

	result += factors.at(0) * results.at(0).at(i) * results.at(0).at(j);
	result += factors.at(1) * results.at(1).at(i) * results.at(1).at(j);
	result += factors.at(2) * results.at(2).at(i) * results.at(2).at(j);
	result += factors.at(3) * results.at(3).at(i) * results.at(3).at(j);

	return result;
}

void test_tensor_decomposition(std::array<T, 4> const& momentum, T denominator)
{
	std::array<std::array<T, 4>, 4> results;
	std::array<T, 4> factors;

	hep::tensor_decompose(momentum, denominator, results, factors);

	// first three vectors must be orthogonal to the momentum
	for (std::size_t i = 0; i != 3; ++i)
	{
		T const product =
			results.at(i)[0] * momentum[0] -
			results.at(i)[1] * momentum[1] -
			results.at(i)[2] * momentum[2] -
			results.at(i)[3] * momentum[3];

		REQUIRE( product < T(1e-9) );
	}

	CHECK_THAT( lhs_matrix_component(momentum, denominator, 0, 0) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 0, 0), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 0, 1) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 0, 1), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 0, 2) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 0, 2), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 0, 3) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 0, 3), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 1, 0) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 1, 0), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 1, 1) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 1, 1), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 1, 2) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 1, 2), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 1, 3) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 1, 3), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 2, 0) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 2, 0), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 2, 1) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 2, 1), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 2, 2) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 2, 2), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 2, 3) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 2, 3), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 3, 0) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 3, 0), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 3, 1) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 3, 1), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 3, 2) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 3, 2), 25) );
	CHECK_THAT( lhs_matrix_component(momentum, denominator, 3, 3) ,
		Catch::WithinULP(rhs_matrix_component(results, factors, 3, 3), 25) );
}

TEST_CASE("test tensor_decompose", "[tensor_decompose]")
{
	// only single space components set, spacelike, denominator larger than one
	test_tensor_decomposition({ T( 0.0), T(1.0), T(0.0), T(0.0) }, T(  5.0));
	test_tensor_decomposition({ T( 0.0), T(0.0), T(1.0), T(0.0) }, T(  5.0));
	test_tensor_decomposition({ T( 0.0), T(0.0), T(0.0), T(1.0) }, T(  5.0));

	// only single space components set, spacelike, denominator smaller than one
	test_tensor_decomposition({ T( 0.0), T(1.0), T(0.0), T(0.0) }, T( -5.0));
	test_tensor_decomposition({ T( 0.0), T(0.0), T(1.0), T(0.0) }, T( -5.0));
	test_tensor_decomposition({ T( 0.0), T(0.0), T(0.0), T(1.0) }, T( -5.0));

	// single space components set, timelike, denominator larger than one
	test_tensor_decomposition({ T( 3.0), T(1.0), T(0.0), T(0.0) }, T(  5.0));
	test_tensor_decomposition({ T( 3.0), T(0.0), T(1.0), T(0.0) }, T(  5.0));
	test_tensor_decomposition({ T( 3.0), T(0.0), T(0.0), T(1.0) }, T(  5.0));

	// single space components set, timelike, denominator larger than one
	test_tensor_decomposition({ T( 3.0), T(1.0), T(0.0), T(0.0) }, T( -5.0));
	test_tensor_decomposition({ T( 3.0), T(0.0), T(1.0), T(0.0) }, T( -5.0));
	test_tensor_decomposition({ T( 3.0), T(0.0), T(0.0), T(1.0) }, T( -5.0));

	// arbitray numbers
	test_tensor_decomposition({ T( 2.0), T(3.0), T(5.0), T(7.0) }, T( 11.0));
	test_tensor_decomposition({ T( 2.0), T(3.0), T(5.0), T(7.0) }, T(-11.0));
	test_tensor_decomposition({ T(-2.0), T(3.0), T(5.0), T(7.0) }, T( 11.0));
	test_tensor_decomposition({ T(-2.0), T(3.0), T(5.0), T(7.0) }, T(-11.0));
}
