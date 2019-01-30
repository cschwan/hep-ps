#include "hep/ps/dipole_vertex.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <cassert>
#include <tuple>

namespace hep
{

bool operator<(dipole_vertex const& a, dipole_vertex const& b)
{
    return std::make_tuple(a.internal(), a.external(), a.unresolved()) <
        std::make_tuple(b.internal(), b.external(), b.unresolved());
}

bool operator==(dipole_vertex const& a, dipole_vertex const& b)
{
    return std::make_tuple(a.internal(), a.external(), a.unresolved()) ==
        std::make_tuple(b.internal(), b.external(), b.unresolved());
}

correction_type correction_type_of(dipole_vertex const& vertex)
{
    if ((vertex.unresolved() == pdg_id_of_photon()) ||
        (vertex.internal() == pdg_id_of_photon()) ||
        (vertex.external() == pdg_id_of_photon()))
    {
        return correction_type::ew;
    }
    else if ((vertex.unresolved() == pdg_id_of_gluon()) ||
        (vertex.internal() == pdg_id_of_gluon()) ||
        (vertex.external() == pdg_id_of_gluon()))
    {
        return correction_type::qcd;
    }

    assert( false );
}

}
