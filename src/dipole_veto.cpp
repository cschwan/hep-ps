#include "hep/ps/dipole_veto.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>

namespace hep
{

dipole_veto default_dipole_veto()
{
	return [](std::vector<int> const& pdg_ids,
		std::vector<final_state> const& final_states, dipole const&)
	{
		auto const& states = pdg_ids_to_states(pdg_ids).second;

		return !std::equal(states.begin(), states.end(), final_states.begin());
	};
}

dipole_veto photon_in_jet()
{
	return [](std::vector<int> const& pdg_ids,
		std::vector<final_state> const& final_states, dipole const&)
	{
		auto const& states = pdg_ids_to_states(pdg_ids).second;
		auto const& mismatch = std::mismatch(final_states.begin(),
			final_states.end(), states.begin());

		if (mismatch.first != final_states.end())
		{
			if ((*mismatch.first != final_state::quark_gluon) ||
				(*mismatch.second != final_state::photon) ||
				!std::equal(mismatch.first + 1, final_states.end(),
				mismatch.second + 1))
			{
				return true;
			}
		}

		// in this case the only difference in the final states is a photon
		// that appears instead of a jet
		return false;
	};
}

}
