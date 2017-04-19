#ifndef HEP_PS_RANDOM_NUMBERS_HPP
#define HEP_PS_RANDOM_NUMBERS_HPP

#include <cassert>
#include <vector>

namespace hep
{

/// A class hold read-only access to a sequence of random numbers.
template <typename T>
class random_numbers
{
public:
	/// Constructor. Note that the given `vector` must be defined in the same
	/// scope or an enclosing one as the object this constructor creates.
	random_numbers(std::vector<T> const& vector)
		: begin(vector.cbegin())
		, end(vector.cend())
	{
	}

	/// Returns the random from the front of the sequence.
	T front()
	{
		assert( begin != end );

		T number = *begin++;

		return number;
	}

	/// Returns the random number from the back of the sequence.
	T back()
	{
		assert( begin != end );

		T number = *--end;

		return number;
	}

private:
	typename std::vector<T>::const_iterator begin;
	typename std::vector<T>::const_iterator end;
};

}

#endif
