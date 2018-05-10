#include "hep/ps/pdg_functions.hpp"

#include <cassert>
#include <stdexcept>
#include <map>

namespace hep
{

final_state pdg_id_to_final_state(int id)
{
	switch (id)
	{
	// quarks
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
	case  -4: return parton::anti_charm;
	case  -3: return parton::anti_strange;
	case  -2: return parton::anti_up;
	case  -1: return parton::anti_down;
	case   1: return parton::down;
	case   2: return parton::up;
	case   3: return parton::strange;
	case   4: return parton::charm;
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
	case parton::anti_charm:   return -4;
	case parton::anti_strange: return -3;
	case parton::anti_up:      return -2;
	case parton::anti_down:    return -1;
	case parton::down:         return  1;
	case parton::up:           return  2;
	case parton::strange:      return  3;
	case parton::charm:        return  4;
	case parton::gluon:        return 21;
	case parton::photon:       return 22;

	default:
		// implementation error
		assert( false );
	}
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

template <typename T>
T pdg_id_to_charge(int id)
{
	switch (id)
	{
	// anti-up-type quarks
	case  -6:
	case  -4:
	case  -2:
		return T(-2.0) / T(3.0);

	// anti-down-type quarks
	case  -5:
	case  -3:
	case  -1:
		return T( 1.0) / T(3.0);

	// down-type quarks
	case   1:
	case   3:
	case   5:
		return T(-1.0) / T(3.0);

	// up-type quarks
	case   2:
	case   4:
	case   6:
		return T( 2.0) / T(3.0);

	// anti-leptons
	case -15:
	case -13:
	case -11:
		return T( 1.0);

	// leptons
	case  11:
	case  13:
	case  15:
		return T(-1.0);

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
	case 22:
		return T();

	default:
		// TODO: NYI
		assert( false );
	}
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template double pdg_id_to_charge(int);

}
