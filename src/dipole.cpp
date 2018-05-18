#include "hep/ps/dipole.hpp"

#include <tuple>

namespace hep
{

bool operator<(dipole const& a, dipole const& b)
{
	return std::make_tuple(
		a.emitter(),
		a.unresolved(),
		a.spectator(),
		a.emitter_type(),
		a.unresolved_type(),
		a.spectator_type(),
		a.corr_type()
	) < std::make_tuple(
		b.emitter(),
		b.unresolved(),
		b.spectator(),
		b.emitter_type(),
		b.unresolved_type(),
		b.spectator_type(),
		b.corr_type()
	);
}

bool operator==(dipole const& a, dipole const& b)
{
	return
		(a.emitter()         == b.emitter())         &&
		(a.unresolved()      == b.unresolved())      &&
		(a.spectator()       == b.spectator())       &&
		(a.emitter_type()    == b.emitter_type())    &&
		(a.unresolved_type() == b.unresolved_type()) &&
		(a.spectator_type()  == b.spectator_type())  &&
		(a.corr_type()       == b.corr_type());
}

}
