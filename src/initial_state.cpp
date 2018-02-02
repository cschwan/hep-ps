#include "hep/ps/initial_state.hpp"

namespace hep
{

parton_set partons_in_initial_state_set(initial_state_set set)
{
	parton_set result;

	for (auto const state : set)
	{
		result.add(state_parton_one(state));
		result.add(state_parton_two(state));
	}

	return result;
}

}
