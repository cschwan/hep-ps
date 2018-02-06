#include "hep/ps/parton.hpp"

#include <cassert>

namespace hep
{

// reference:
// http://pdg.lbl.gov/2017/reviews/rpp2017-rev-monte-carlo-numbering.pdf
int parton_to_pdg_id(parton p)
{
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
