#include "hep/ps/vertex.hpp"

#include <algorithm>

namespace hep
{

vertex::vertex(
	std::size_t p1,
	std::size_t p2,
	std::size_t p3,
	std::size_t p4
)
	: p_{p1, p2, p3, p4}
{
}

std::array<std::size_t, 4> const& vertex::p() const
{
	return p_;
}

}
