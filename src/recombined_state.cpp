#include "hep/ps/recombined_state.hpp"

#include <cassert>

namespace hep
{

recombined_state trivially_recombine(final_state state)
{
    switch (state)
    {
    case final_state::charged_lepton:
        return recombined_state::dressed_lepton;

    case final_state::quark_gluon:
        return recombined_state::jet;

    case final_state::photon:
        return recombined_state::isolated_photon;

    case final_state::neutrino:
        return recombined_state::missing_momentum;

    default:
        // this means I didn't cover all the cases
        assert( false );
    }
}

}
