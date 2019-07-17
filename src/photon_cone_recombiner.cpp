#include "hep/ps/photon_cone_recombiner.hpp"
#include "hep/ps/psp.hpp"
#include "hep/ps/psp_type.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <limits>

namespace hep
{

template <typename T>
photon_cone_recombiner<T>::photon_cone_recombiner(T radius)
    : radius2_{radius * radius}
{
}

template <typename T>
void photon_cone_recombiner<T>::recombine(
    std::vector<T> const& phase_space,
    std::vector<final_state> const& final_states,
    std::vector<T>& recombined_phase_space,
    std::vector<recombined_state>& recombined_states
) {
    candidates_.clear();
    recombined_states.clear();

    std::size_t j;
    bool found_photon = false;

    for (std::size_t i = 0; i != final_states.size(); ++i)
    {
        switch (final_states.at(i))
        {
        case final_state::charged_lepton:
            recombined_states.push_back(recombined_state::dressed_lepton);
            break;

        case final_state::quark_gluon:
            recombined_states.push_back(recombined_state::jet);
            candidates_.push_back(i + 2);
            break;

        case final_state::photon:
            recombined_states.push_back(recombined_state::isolated_photon);
            // TODO: implement general case
            assert( !found_photon );
            found_photon = true;
            j = i + 2;
            break;

        case final_state::neutrino:
            recombined_states.push_back(recombined_state::missing_momentum);
            break;
        }
    }

    recombined_phase_space = phase_space;

    if (!found_photon)
    {
        return;
    }

    std::size_t min_i = 0;
    T min_radius2 = std::numeric_limits<T>::max();

    psp<T> ps{recombined_phase_space, recombined_states, T(), psp_type::pos_rap};

    for (std::size_t i = 0; i != final_states.size(); ++i)
    {
        if ((final_states.at(i) == final_state::photon) ||
            (final_states.at(i) == final_state::neutrino))
        {
            continue;
        }

        T r2ij = ps.dist2(i + 2, j);

        if (r2ij < min_radius2)
        {
            min_radius2 = r2ij;
            min_i = i + 2;
        }
    }

    if (min_radius2 < radius2_)
    {
        recombined_states.erase(std::next(recombined_states.begin(), j - 2));

        recombined_phase_space.at(4 * min_i + 0) += recombined_phase_space.at(4 * j + 0);
        recombined_phase_space.at(4 * min_i + 1) += recombined_phase_space.at(4 * j + 1);
        recombined_phase_space.at(4 * min_i + 2) += recombined_phase_space.at(4 * j + 2);
        recombined_phase_space.at(4 * min_i + 3) += recombined_phase_space.at(4 * j + 3);

        recombined_phase_space.erase(std::next(recombined_phase_space.begin(), 4 * j),
            std::next(recombined_phase_space.begin(), 4 * j + 4));
    }
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class photon_cone_recombiner<double>;

}
