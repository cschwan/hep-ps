#include "hep/ps/dipole.hpp"

#include <cassert>
#include <tuple>

namespace hep
{

dipole::dipole(
	std::size_t emitter,
	std::size_t unresolved,
	std::size_t spectator,
	particle_type emitter_type,
	particle_type unresolved_type,
	particle_type spectator_type,
	correction_type corr_type
)
	: emitter_(emitter)
	, unresolved_(unresolved)
	, spectator_(spectator)
	, emitter_type_(emitter_type)
	, unresolved_type_(unresolved_type)
	, spectator_type_(spectator_type)
	, corr_type_(corr_type)
{
	int type = (emitter < 2) | ((spectator < 2) << 1);

	switch (type)
	{
	case 0:
		type_ = dipole_type::final_final;
		break;

	case 1:
		type_ = dipole_type::initial_final;
		break;

	case 2:
		type_ = dipole_type::final_initial;
		break;

	case 3:
		type_ = dipole_type::initial_initial;
		break;

	default:
		assert( false );
	}
}

bool operator<(dipole const& a, dipole const& b)
{
	return std::make_tuple(
		a.corr_type(),
		a.emitter(),
		a.unresolved(),
		a.spectator(),
		a.emitter_type(),
		a.unresolved_type(),
		a.spectator_type()
	) < std::make_tuple(
		b.corr_type(),
		b.emitter(),
		b.unresolved(),
		b.spectator(),
		b.emitter_type(),
		b.unresolved_type(),
		b.spectator_type()
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

namespace std
{

size_t hash<hep::dipole>::operator()(hep::dipole const& dipole) const
{
	// we assume that most dipoles don't have indices larger than 15 (five
	// bits) and are basically defined by those three indices
	return
		((dipole.emitter()    & 0b11111) <<  0) ||
		((dipole.spectator()  & 0b11111) <<  5) ||
		((dipole.unresolved() & 0b11111) << 10);
}

}
