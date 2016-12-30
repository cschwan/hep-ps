#include "hep/ps/fortran_helper.hpp"

#include <cassert>

namespace hep
{

template <typename T>
void fortran_ordering_to_cpp(std::vector<T>& momenta)
{
	// size must be a multiple of four
	assert( momenta.size() % 4 == 0 );

	std::size_t const n = momenta.size() / 4;
	std::vector<T> buffer = momenta;

	for (std::size_t i = 0; i != buffer.size(); ++i)
	{
		// component index
		std::size_t const k = i / n;
		// particle index
		std::size_t const j = i % n;

		momenta.at(4 * j + k) = buffer.at(i);
	}
}

template <typename T>
void cpp_ordering_to_fortran(std::vector<T>& momenta)
{
	// size must be a multiple of four
	assert( momenta.size() % 4 == 0 );

	std::size_t const n = momenta.size() / 4;
	std::vector<T> buffer = momenta;

	for (std::size_t i = 0; i != buffer.size(); ++i)
	{
		// component index
		std::size_t const k = i % 4;
		// particle index
		std::size_t const j = i / 4;

		momenta.at(n * k + j) = buffer.at(i);
	}
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template void fortran_ordering_to_cpp(std::vector<double>&);
template void cpp_ordering_to_fortran(std::vector<double>&);

}
