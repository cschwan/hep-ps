#include "hep/ps/int_dipole.hpp"

namespace hep
{

bool operator<(int_dipole const& a, int_dipole const& b)
{
    if (a.type() == insertion_term_type::born)
    {
        if (b.type() == insertion_term_type::born)
        {
            return a.initial_particle() < b.initial_particle();
        }

        return true;
    }
    else if (b.type() == insertion_term_type::born)
    {
        return false;
    }

    return std::make_tuple(a.emitter(), a.spectator(), a.vertex()) <
        std::make_tuple(b.emitter(), b.spectator(), b.vertex());
}

bool operator==(int_dipole const& a, int_dipole const& b)
{
    return (a.vertex() == b.vertex()) && (a.type() == b.type()) &&
        ((a.type() == insertion_term_type::born) ? (a.initial_particle() == b.initial_particle())
        : ((a.emitter() == b.emitter()) && (a.spectator() == b.spectator())));
}


}

namespace std
{

size_t hash<hep::int_dipole>::operator()(hep::int_dipole const& dip) const
{
    // TODO: take `vertex()` into account

    // we assume that most dipoles don't have indices larger than 15 (five
    // bits) and are basically defined by those three indices
    return (dip.type() == hep::insertion_term_type::born) ?
        ((dip.initial_particle() & 0b11111) <<  0) :
        (((dip.emitter()         & 0b11111) <<  0) ||
         ((dip.spectator()       & 0b11111) <<  5));
}

}
