#include "hep/ps/permutation.hpp"

#include "catch2/catch.hpp"

#include <array>

TEST_CASE("check inverse permutation", "[permutation]")
{
    constexpr int s  = 0;
    constexpr int u  = 1;
    constexpr int ve = 2;
    constexpr int el = 3;
    constexpr int vm = 4;
    constexpr int mu = 5;
    constexpr int d  = 6;
    constexpr int c  = 7;
    constexpr int g  = 8;

    std::array<int, 9> permutation;
    std::array<int, 9> inverse;

    permutation = { s, g, ve, el, vm, mu, d, c, u };
    inverse = hep::inverse_permutation(permutation);

    // `permutation` is self-inverse
    CHECK( inverse[0] == permutation[0] );
    CHECK( inverse[1] == permutation[1] );
    CHECK( inverse[2] == permutation[2] );
    CHECK( inverse[3] == permutation[3] );
    CHECK( inverse[4] == permutation[4] );
    CHECK( inverse[5] == permutation[5] );
    CHECK( inverse[6] == permutation[6] );
    CHECK( inverse[7] == permutation[7] );
    CHECK( inverse[8] == permutation[8] );

    permutation = { u, c, ve, el, vm, mu, d, s, g };
    inverse = hep::inverse_permutation(permutation);

    CHECK( inverse[0] == 7 );
    CHECK( inverse[1] == 0 );
    CHECK( inverse[2] == 2 );
    CHECK( inverse[3] == 3 );
    CHECK( inverse[4] == 4 );
    CHECK( inverse[5] == 5 );
    CHECK( inverse[6] == 6 );
    CHECK( inverse[7] == 1 );
    CHECK( inverse[8] == 8 );
}
