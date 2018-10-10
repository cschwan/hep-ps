#include "hep/ps/trivial_recombiner.hpp"

#include <algorithm>

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
    recombined_states.resize(final_states.size());

    std::transform(final_states.begin(), final_states.end(), recombined_states.begin(),
        trivially_recombine);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class trivial_recombiner<double>;

}
