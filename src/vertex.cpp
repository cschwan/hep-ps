#include "hep/ps/vertex.hpp"

#include <algorithm>

namespace hep
{

vertex::vertex(
	inv_idx const& inv1,
	inv_idx const& inv2,
	inv_idx const& inv3,
	inv_idx const& inv4
)
	: invariants_{inv1, inv2, inv3, inv4}
{
}

std::array<inv_idx, 4> const& vertex::invariants() const
{
	return invariants_;
}

}
