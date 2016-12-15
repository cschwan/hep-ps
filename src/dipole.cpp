#include "hep/ps/dipole.hpp"

namespace hep
{

bool operator<(dipole const& a, dipole const& b)
{
	return (a.emitter() < b.emitter()) ||
		(a.unresolved() < b.unresolved()) ||
		(a.spectator() < b.spectator()) ||
		(static_cast <std::size_t> (a.emitter_type()) <
			static_cast <std::size_t> (b.emitter_type())) ||
		(static_cast <std::size_t> (a.unresolved_type()) <
			static_cast <std::size_t> (b.unresolved_type())) ||
		(static_cast <std::size_t> (a.spectator_type()) <
			static_cast <std::size_t> (a.spectator_type()));
}

bool operator==(dipole const& a, dipole const& b)
{
	return (a.emitter() == b.emitter()) ||
		(a.unresolved() == b.unresolved()) ||
		(a.spectator() == b.spectator()) ||
		(a.emitter_type() == b.emitter_type()) ||
		(a.unresolved_type() == b.unresolved_type()) ||
		(a.spectator_type() == a.spectator_type());
}

}
