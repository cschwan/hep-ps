#include "hep/ps/diagram.hpp"

namespace hep
{

diagram::diagram(
	std::vector<vertex> const& vertices,
	std::vector<propagator> const& propagators
)
	: vertices_{vertices}
	, propagators_{propagators}
{
}

std::vector<vertex> const& diagram::vertices() const
{
	return vertices_;
}

std::vector<propagator> const& diagram::propagators() const
{
	return propagators_;
}

}
