#include "hep/ps/pdg_functions.hpp"

#include <cassert>

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

initial_state pdg_ids_to_initial_state(int id1, int id2)
{
	if (id1 == 2 && id2 == 2)
	{
		return initial_state::q43_uu;
	}
	else if (id1 == -3 && id2 == 2)
	{
		return initial_state::q33_su;
	}
	else
	{
		// TODO: NYI
		assert( false );
	}
}

}
