#include "hep/ps/pdg_functions.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <map>

namespace hep
{

final_state pdg_id_to_final_state(int id)
{
	switch (id)
	{
	case  -6:
	case  -5:
	case  -4:
	case  -3:
	case  -2:
	case  -1:
	case   1:
	case   2:
	case   3:
	case   4:
	case   5:
	case   6:
	case  21:
		return final_state::quark_gluon;

	case -15:
	case -13:
	case -11:
	case  11:
	case  13:
	case  15:
		return final_state::charged_lepton;

	case -16:
	case -14:
	case -12:
	case  12:
	case  14:
	case  16:
		return final_state::neutrino;

	case 22:
		return final_state::photon;

	default:
		// TODO: NYI
		assert( false );
	}
}

parton pdg_id_to_parton(int id)
{
	switch (id)
	{
	case  -5: return parton::anti_bottom;
	case  -4: return parton::anti_charm;
	case  -3: return parton::anti_strange;
	case  -2: return parton::anti_up;
	case  -1: return parton::anti_down;
	case   1: return parton::down;
	case   2: return parton::up;
	case   3: return parton::strange;
	case   4: return parton::charm;
	case   5: return parton::bottom;
	case  21: return parton::gluon;
	case  22: return parton::photon;

	default:
		// TODO: NYI
		assert( false );
	}
}

particle_type pdg_id_to_particle_type(int id)
{
	switch (id)
	{
	case -16:
	case -15:
	case -14:
	case -13:
	case -12:
	case -11:
	case  -6:
	case  -5:
	case  -4:
	case  -3:
	case  -2:
	case  -1:
	case   1:
	case   2:
	case   3:
	case   4:
	case   5:
	case   6:
	case  11:
	case  12:
	case  13:
	case  14:
	case  15:
	case  16:
		return particle_type::fermion;

	case  22:
	case  21:
		return particle_type::boson;

	default:
		// TODO: NYI
		assert( false );
	}
}

int parton_to_pdg_id(parton p)
{
	// reference:
	// http://pdg.lbl.gov/2017/reviews/rpp2017-rev-monte-carlo-numbering.pdf
	switch (p)
	{
	case parton::anti_bottom:  return -5;
	case parton::anti_charm:   return -4;
	case parton::anti_strange: return -3;
	case parton::anti_up:      return -2;
	case parton::anti_down:    return -1;
	case parton::down:         return  1;
	case parton::up:           return  2;
	case parton::strange:      return  3;
	case parton::charm:        return  4;
	case parton::bottom:       return  5;
	case parton::gluon:        return 21;
	case parton::photon:       return 22;

	default:
		// implementation error
		assert( false );
	}
}

bool pdg_id_has_color(int id)
{
	switch (id)
	{
	case  -6:
	case  -5:
	case  -4:
	case  -3:
	case  -2:
	case  -1:
	case   1:
	case   2:
	case   3:
	case   4:
	case   5:
	case   6:
	case  21:
		return true;

	case -16:
	case -15:
	case -14:
	case -13:
	case -12:
	case -11:
	case  11:
	case  13:
	case  15:
	case  12:
	case  14:
	case  16:
	case  22:
		return false;

	default:
		// TODO: NYI
		assert( false );
	}
}

bool pdg_id_has_charge(int id)
{
	switch (id)
	{
	case -15:
	case -13:
	case -11:
	case  -6:
	case  -5:
	case  -4:
	case  -3:
	case  -2:
	case  -1:
	case   1:
	case   2:
	case   3:
	case   4:
	case   5:
	case   6:
	case  11:
	case  13:
	case  15:
		return true;

	// neutrinos
	case -16:
	case -14:
	case -12:
	case  12:
	case  14:
	case  16:
	case  21:
	case  22:
		return false;

	default:
		// TODO: NYI
		assert( false );
	}
}

int pdg_id_of_photon()
{
	return 22;
}

int pdg_id_of_gluon()
{
	return 21;
}

std::vector<int> ol_process_string_to_pdg_ids(std::string const& process)
{
	// format int int `->` int int ...
	std::vector<int> result;

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

		result.push_back(std::stoi(substring));
	}

	return result;
}

std::string pdg_ids_to_ol_process_string(std::vector<int> const& ids)
{
	std::string result;

	// reserve space for each id and the separator
	result.reserve((ids.size() + 1) * 4);

	result.append(std::to_string(ids.at(0)));
	result.append(" ");
	result.append(std::to_string(ids.at(1)));
	result.append(" ->");

	for (std::size_t i = 2; i != ids.size(); ++i)
	{
		result.append(" ");
		result.append(std::to_string(ids.at(i)));
	}

	return result;
}

bool pdg_id_is_gluon(int id)
{
	return id == pdg_id_of_gluon();
}

bool pdg_id_is_quark(int id)
{
	switch (id)
	{
	case -6:
	case -5:
	case -4:
	case -3:
	case -2:
	case -1:
	case  1:
	case  2:
	case  3:
	case  4:
	case  5:
	case  6:
		return true;

	default:
		return false;
	}
}

bool pdg_id_is_photon(int id)
{
	return id == pdg_id_of_photon();
}

std::size_t final_state_symmetry_factor(std::vector<int> const& ids)
{
	std::size_t result = 1;
	std::vector<int> copy(ids.begin() + 2, ids.end());

	while (!copy.empty())
	{
		int id = copy.front();
		std::size_t count = std::count(copy.begin(), copy.end(), id);
		copy.erase(std::remove(copy.begin(), copy.end(), id), copy.end());

		while (count > 1)
		{
			result *= count;
			--count;
		}
	}

	return result;
}

int pdg_id_to_charge_times_three(int id)
{
	switch (id)
	{
	// anti-up-type quarks
	case  -6:
	case  -4:
	case  -2:
		return -2;

	// anti-down-type quarks
	case  -5:
	case  -3:
	case  -1:
		return 1;

	// down-type quarks
	case   1:
	case   3:
	case   5:
		return -1;

	// up-type quarks
	case   2:
	case   4:
	case   6:
		return 2;

	// anti-leptons
	case -15:
	case -13:
	case -11:
		return 3;

	// leptons
	case  11:
	case  13:
	case  15:
		return -3;

	// neutrinos
	case -16:
	case -14:
	case -12:
	case  12:
	case  14:
	case  16:
	// gluon
	case  21:
	// photon
	case  22:
		return 0;

	default:
		// TODO: NYI
		assert( false );
	}
}

std::pair<initial_state, std::vector<final_state>> pdg_ids_to_states(
	std::vector<int> const& ids
) {
	std::vector<final_state> final_states(ids.size() - 2);
	std::transform(ids.begin() + 2, ids.end(), final_states.begin(),
		pdg_id_to_final_state);

	auto const state = partons_to_initial_state(pdg_id_to_parton(ids.at(0)),
		pdg_id_to_parton(ids.at(1)));

	return { state , final_states };
}

}
