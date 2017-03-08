#ifndef HEP_PS_INITIAL_STATE_SET_HPP
#define HEP_PS_INITIAL_STATE_SET_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016-2017  Christopher Schwan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hep/ps/initial_state.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <type_traits>
#include <utility>

namespace hep
{

class initial_state_set
{
public:
	class const_iterator
	{
	public:
		using iterator_category = std::input_iterator_tag;
		using value_type = initial_state;
		using difference_type = std::ptrdiff_t;
		using pointer = initial_state const*;
		using reference = initial_state;

		const_iterator()
			: set_{}
		{
		}

		const_iterator(const_iterator const&) = default;

		~const_iterator() = default;

		const_iterator& operator=(const const_iterator&) = default;

		explicit const_iterator(std::size_t set)
			: set_(set)
		{
		}

		// Swappable
		friend void swap(const_iterator& a, const_iterator& b)
		{
			std::swap(a.set_, b.set_);
		}

		// Iterator
		const_iterator& operator++()
		{
			std::size_t const least_significant_bit = set_ &
				-static_cast <std::make_signed<std::size_t>::type> (set_);
			set_ ^= least_significant_bit;

			return *this;
		}

		const_iterator operator++(int)
		{
			auto const result = *this;
			++(*this);
			return result;
		}

		bool operator==(const_iterator other) const
		{
			return set_ == other.set_;
		}

		bool operator!=(const_iterator other) const
		{
			return !(*this == other);
		}

		reference operator*() const
		{
			// the following algorithm works only for 31-bit integers (the 32nd
			// bit is needed for the sign)
			assert( set_ <= std::numeric_limits<std::int32_t>::max() );

			// algorithm copied from "Bit Twiddling Hacks"
			std::size_t v = set_;
			std::size_t c = 32;
			v &= -static_cast <std::make_signed<std::size_t>::type> (v);

			if (v) c--;
			if (v & 0x0000FFFF) c -= 16;
			if (v & 0x00FF00FF) c -= 8;
			if (v & 0x0F0F0F0F) c -= 4;
			if (v & 0x33333333) c -= 2;
			if (v & 0x55555555) c -= 1;

			return static_cast <initial_state> (c);
		}

	private:
		std::size_t set_;
	};

	const_iterator begin() const
	{
		return const_iterator(set_);
	}

	const_iterator end() const
	{
		return const_iterator();
	}

	initial_state_set()
		: set_{}
	{
	}

	initial_state_set(std::initializer_list<initial_state> list)
		: set_(std::accumulate(list.begin(), list.end(), 0u,
			[](std::size_t set, initial_state state) {
				return set | (1 << static_cast <std::size_t> (state));
		  }))
	{
	}

	bool empty() const
	{
		return set_ == 0;
	}

	void subtract(initial_state_set other_set)
	{
		set_ ^= (set_ & other_set.set_);
	}

	void add(initial_state state)
	{
		auto const index = static_cast <std::size_t> (state);
		set_ |= (1 << index);
	}

	bool includes(initial_state state) const
	{
		auto const index = static_cast <std::size_t> (state);
		return set_ & (1 << index);
	}

	bool operator==(initial_state_set other) const
	{
		return set_ == other.set_;
	}

private:
	std::size_t set_;
};

}

#endif
