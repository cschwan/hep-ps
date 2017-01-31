#ifndef HEP_PS_PERMUTATION_HPP
#define HEP_PS_PERMUTATION_HPP

#include <array>
#include <cstddef>

namespace hep
{

template <typename T, std::size_t N>
inline std::array<T, N> inverse_permutation(std::array<T, N> permutation)
{
	// TODO: check if `permutation` is actually a permutation?

	std::array<T, N> inverse;

	for (std::size_t i = 0; i != N; ++i)
	{
		inverse.at(permutation[i]) = i;
	}

	return inverse;
}

}

#endif
