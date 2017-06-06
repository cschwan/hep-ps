#include "hep/ps/propagator.hpp"

namespace hep
{

propagator::propagator(std::size_t type, inv_idx const& invariant)
	: type_{type}
	, invariant_{invariant}
{
}

std::size_t propagator::type() const
{
	return type_;
}

inv_idx const& propagator::invariant() const
{
	return invariant_;
}


}
