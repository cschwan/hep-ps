#include "hep/ps/initial_state.hpp"

#include <map>
#include <stdexcept>
#include <utility>

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

initial_state partons_to_initial_state(parton one, parton two)
{
	static std::map<std::pair<parton, parton>, initial_state> map;
	static bool initialized = false;

	if (!initialized)
	{
		for (auto state : initial_state_list())
		{
			map.emplace(std::make_pair(state_parton_one(state),
				state_parton_two(state)), state);
		}

		initialized = true;
	}

	auto const result = map.find({ one, two });

	if (result == map.end())
	{
		throw std::invalid_argument("no initial state found for given partons");
	}

	return result->second;
}

template <typename T>
T casimir_operator(initial_state state, std::size_t index)
{
	if (index == 0)
	{
		switch (state_parton_one(state))
		{
		case parton::photon:  return T();
		case parton::gluon:   return T(3.0);
		default:              return T(4.0) / T(3.0);
		}
	}
	else if (index == 1)
	{
		switch (state_parton_two(state))
		{
		case parton::photon:  return T();
		case parton::gluon:   return T(3.0);
		default:              return T(4.0) / T(3.0);
		}
	}

	assert( false );
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template double casimir_operator(initial_state, std::size_t);

}
