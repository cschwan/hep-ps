#include "hep/ps/silence_third_party_libs.hpp"

#include <LHAPDF/PDF.h>

namespace hep
{

void silence_third_party_libs()
{
	LHAPDF::setVerbosity(0);
}

}
