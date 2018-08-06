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
	else if (b.type() == insertion_term_type::born)
	{
		return false;
	}

	return std::make_tuple(
		a.emitter(),
		a.emitter_type(),
		a.spectator()
	) < std::make_tuple(
		b.emitter(),
		b.emitter_type(),
		b.spectator()
	);
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

namespace std
{

size_t hash<hep::insertion_term>::operator()(
	hep::insertion_term const& term
) const {
	// we assume that most dipoles don't have indices larger than 15 (five
	// bits) and are basically defined by those three indices
	return (term.type() == hep::insertion_term_type::born) ?
		((term.initial_particle() & 0b11111) <<  0) :
		(((term.emitter()         & 0b11111) <<  0) ||
		 ((term.spectator()       & 0b11111) <<  5));
}

}
