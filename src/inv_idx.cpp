#include "hep/ps/inv_idx.hpp"

namespace hep
{

inv_idx::inv_idx()
	: representation_{0}
	, n_{0}
{
}

inv_idx::inv_idx(std::size_t particle_index, std::size_t n)
	: representation_{particle_index}
	, n_{n}
{
}

std::size_t inv_idx::representation() const
{
	return representation_;
}

std::size_t inv_idx::n() const
{
	return n_;
}

bool is_s_like(inv_idx const& index)
{
	bool const has_one = index.representation() & 1;
	bool const has_two = index.representation() & 2;

	// an s-like invariant must have either both or none of the initial state
	// momenta
	return has_one == has_two;
}

bool is_t_like(inv_idx const& index)
{
	return !is_s_like(index);
}

}
