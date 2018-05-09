#include "hep/ps/ol_helpers.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>

namespace hep
{

std::pair<initial_state, std::vector<final_state>> parse_ol_process_string(
	std::string const& process
) {
	// format int int `->` int int ...
	std::vector<int> pdg_ids;

	std::size_t begin = 0;
	std::size_t end = 0;

	while (end != process.size())
	{
		begin = process.find_first_not_of(' ', begin);
		end = process.find(' ', begin + 1);

		if (end == std::string::npos)
		{
			end = process.size();
		}

		auto const& substring = process.substr(begin, end - begin);

		begin = end;

		if (substring == "->")
		{
			continue;
		}

		pdg_ids.push_back(std::stoi(substring));
	}

	std::vector<final_state> final_states(pdg_ids.size() - 2);
	std::transform(pdg_ids.begin() + 2, pdg_ids.end(), final_states.begin(),
		[](int id) { return pdg_id_to_final_state(id); });

	return { pdg_ids_to_initial_state(pdg_ids.at(0), pdg_ids.at(1)),
		final_states };
}

}
