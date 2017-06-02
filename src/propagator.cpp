#include "hep/ps/propagator.hpp"

namespace hep
{

propagator::propagator(std::size_t type, std::size_t momenta)
	: type_{type}
	, momenta_{momenta}
{
}

std::size_t propagator::type() const
{
	return type_;
}

std::size_t propagator::momenta() const
{
	return momenta_;
}


}
