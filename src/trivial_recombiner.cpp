#include "hep/ps/trivial_recombiner.hpp"

namespace hep
{

template <typename T>
void trivial_recombiner<T>::recombine(
    std::vector<T> const& phase_space,
    std::vector<final_state> const& final_states,
    std::vector<T>& recombined_phase_space,
    std::vector<recombined_state>& recombined_states
) const {
    recombined_phase_space = phase_space;
    recombined_states.clear();

    for (auto const state : final_states)
    {
        switch (state)
        {
        case final_state::charged_lepton:
            recombined_states.push_back(recombined_state::dressed_lepton);
            break;

        case final_state::quark_gluon:
            recombined_states.push_back(recombined_state::jet);
            break;

        case final_state::photon:
            recombined_states.push_back(recombined_state::isolated_photon);
            break;

        case final_state::neutrino:
            recombined_states.push_back(recombined_state::missing_momentum);
            break;
        }
    }
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class trivial_recombiner<double>;

}
