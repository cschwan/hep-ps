#include "hep/ps/suppress_banners.hpp"

namespace hep
{

void suppress_banners(bool suppress)
{
    suppress_banners() = suppress;
}

bool& suppress_banners()
{
    static bool suppression = false;
    return suppression;
}

}
