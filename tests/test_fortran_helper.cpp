#include "hep/ps/fortran_helper.hpp"

#include <catch.hpp>

#include <vector>

using T = HEP_TYPE_T;

TEST_CASE("FORTRAN -> C++", "[fortran_helper]")
{
    // buffer is in FORTRAN ordering
    std::vector<T> buffer{
        // 1        2        3
        T( 1.0), T( 2.0), T( 3.0),
        T( 4.0), T( 5.0), T( 6.0),
        T( 7.0), T( 8.0), T( 9.0),
        T(10.0), T(11.0), T(12.0)
    };

    hep::fortran_ordering_to_cpp(buffer);

    // first momentum
    REQUIRE( buffer.at( 0) == T( 1.0) );
    REQUIRE( buffer.at( 1) == T( 4.0) );
    REQUIRE( buffer.at( 2) == T( 7.0) );
    REQUIRE( buffer.at( 3) == T(10.0) );
    // second momentum
    REQUIRE( buffer.at( 4) == T( 2.0) );
    REQUIRE( buffer.at( 5) == T( 5.0) );
    REQUIRE( buffer.at( 6) == T( 8.0) );
    REQUIRE( buffer.at( 7) == T(11.0) );
    // third momentum
    REQUIRE( buffer.at( 8) == T( 3.0) );
    REQUIRE( buffer.at( 9) == T( 6.0) );
    REQUIRE( buffer.at(10) == T( 9.0) );
    REQUIRE( buffer.at(11) == T(12.0) );
}

TEST_CASE("C++ -> FORTRAN", "[fortran helper]")
{
    // buffer is in C++ ordering
    std::vector<T> buffer{
        /* 1 */ T( 1.0), T( 2.0), T( 3.0), T( 4.0),
        /* 2 */ T( 5.0), T( 6.0), T( 7.0), T( 8.0),
        /* 3 */ T( 9.0), T(10.0), T(11.0), T(12.0)
    };

    hep::cpp_ordering_to_fortran(buffer);

    // energies
    REQUIRE( buffer.at( 0) == T( 1.0) );
    REQUIRE( buffer.at( 1) == T( 5.0) );
    REQUIRE( buffer.at( 2) == T( 9.0) );
    // x-components
    REQUIRE( buffer.at( 3) == T( 2.0) );
    REQUIRE( buffer.at( 4) == T( 6.0) );
    REQUIRE( buffer.at( 5) == T(10.0) );
    // y-components
    REQUIRE( buffer.at( 6) == T( 3.0) );
    REQUIRE( buffer.at( 7) == T( 7.0) );
    REQUIRE( buffer.at( 8) == T(11.0) );
    // z-components
    REQUIRE( buffer.at( 9) == T( 4.0) );
    REQUIRE( buffer.at(10) == T( 8.0) );
    REQUIRE( buffer.at(11) == T(12.0) );
}
