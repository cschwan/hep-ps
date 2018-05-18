#include "hep/ps/insertion_term.hpp"

namespace hep
{

bool operator<(insertion_term const& a, insertion_term const& b)
{
	if (a.corr_type() != b.corr_type())
	{
		return static_cast <int> (a.corr_type()) <
			static_cast <int> (b.corr_type());
	}

	if (a.type() == insertion_term_type::born)
	{
		if (b.type() == insertion_term_type::born)
		{
			return a.initial_particle() < b.initial_particle();
		}

		return true;
	}
	else if (b.type() != insertion_term_type::born)
	{
		return (a.emitter() < b.emitter()) ||
			(static_cast <int> (a.emitter_type()) <
				static_cast <int> (b.emitter_type())) ||
			(a.spectator() < b.spectator());
	}

	return false;
}

bool operator==(insertion_term const& a, insertion_term const& b)
{
	return (a.corr_type() == b.corr_type()) && (a.type() == b.type()) &&
		((a.type() == insertion_term_type::born)
		? (a.initial_particle() == b.initial_particle())
		: ((a.emitter() == b.emitter()) &&
			(a.emitter_type() == b.emitter_type()) &&
			(a.spectator() == b.spectator())));
}

}
