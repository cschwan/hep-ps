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

}
