#ifndef HEP_PS_ENUM_HPP
#define HEP_PS_ENUM_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017  Christopher Schwan
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

#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

/// Returns the number of arguments given to this preprocessor macro.
#define HEP_ENUM_SIZEOF(...) \
	std::tuple_size<decltype (std::make_tuple(__VA_ARGS__))>::value

// TODO: make the `enum` an `enum class`

/// Defines an enumeration `name` and a `constexpr` function `name_list()`. The
/// arguments for this macro after `name` are the possible values of the
/// enumeration. The function returns a `std::array<name, ...>` with the
/// possible values of the enumeration.
#define HEP_ENUM(name, ...)                                                    \
	enum name : std::size_t                                                    \
	{                                                                          \
		__VA_ARGS__                                                            \
	};                                                                         \
	constexpr std::array<name, HEP_ENUM_SIZEOF(__VA_ARGS__)> name ## _list()   \
	{                                                                          \
		return { __VA_ARGS__ };                                                \
	}                                                                          \
	enum name : std::size_t

/// Defines a template class `name_array` for a previously defined enumeration
/// defined with \ref HEP_ENUM. The class `name_array<T>` is a wrapper over
/// `std::array<T, ...>` that supports access to its values by keys that are
/// values of the enumeration. A default constructor is provided that
/// zero-initializes all elements of the array.
#define HEP_ENUM_ARRAY(name)                                                   \
	template <typename T>                                                      \
	class name ## _array                                                       \
	{                                                                          \
	public:                                                                    \
		/** Default constructor */                                             \
		name ## _array()                                                       \
			: array_{{}}                                                       \
		{                                                                      \
		}                                                                      \
		                                                                       \
		/** Uses `index` to access the corresponding element. */               \
		T& operator[](name index)                                              \
		{                                                                      \
			return array_[static_cast <std::size_t> (index)];                  \
		}                                                                      \
		                                                                       \
		/** Uses `index` to access the corresponding element. */               \
		T operator[](name index) const                                         \
		{                                                                      \
			return array_[static_cast <std::size_t> (index)];                  \
		}                                                                      \
		                                                                       \
		/** Adds the elements of `other` this objects' elements. */            \
		name ## _array& operator+=(name ## _array const& other)                \
		{                                                                      \
			for (std::size_t i = 0; i != array_.size(); ++i)                   \
			{                                                                  \
				array_[i] += other.array_[i];                                  \
			}                                                                  \
	                                                                           \
			return *this;                                                      \
		}                                                                      \
	                                                                           \
	private:                                                                   \
		std::array<T, name ## _list().size()> array_;                          \
	}

///
#define HEP_ENUM_MAP(name)                                                     \
	template <typename T>                                                      \
	using name ## _map = std::vector<std::pair<name, T>>

/// Defines a class that is able to contain a unique element of the given
/// enumeration `name` that must be previously declared with \ref HEP_ENUM.
#define HEP_ENUM_SET(name)                                                     \
	class name ## _set                                                         \
	{                                                                          \
	public:                                                                    \
		class const_iterator                                                   \
		{                                                                      \
		public:                                                                \
			using iterator_category = std::input_iterator_tag;                 \
			using value_type        = name;                                    \
			using difference_type   = std::ptrdiff_t;                          \
			using pointer           = name const*;                             \
			using reference         = name;                                    \
                                                                               \
			const_iterator()                                                   \
				: set_{}                                                       \
			{                                                                  \
			}                                                                  \
                                                                               \
			const_iterator(const_iterator const&) = default;                   \
                                                                               \
			~const_iterator() = default;                                       \
                                                                               \
			const_iterator& operator=(const const_iterator&) = default;        \
                                                                               \
			explicit const_iterator(std::size_t set)                           \
				: set_{set}                                                    \
			{                                                                  \
			}                                                                  \
                                                                               \
			friend void swap(const_iterator& a, const_iterator& b)             \
			{                                                                  \
				std::swap(a.set_, b.set_);                                     \
			}                                                                  \
                                                                               \
			const_iterator& operator++()                                       \
			{                                                                  \
				std::size_t const least_significant_bit = set_ &               \
					-static_cast <std::make_signed<std::size_t>::type> (set_); \
				set_ ^= least_significant_bit;                                 \
                                                                               \
				return *this;                                                  \
			}                                                                  \
                                                                               \
			const_iterator operator++(int)                                     \
			{                                                                  \
				auto const result = *this;                                     \
				++(*this);                                                     \
				return result;                                                 \
			}                                                                  \
                                                                               \
			bool operator==(const_iterator other) const                        \
			{                                                                  \
				return set_ == other.set_;                                     \
			}                                                                  \
                                                                               \
			bool operator!=(const_iterator other) const                        \
			{                                                                  \
				return !(*this == other);                                      \
			}                                                                  \
                                                                               \
			reference operator*() const                                        \
			{                                                                  \
				assert( set_ <= std::numeric_limits<std::int32_t>::max() );    \
                                                                               \
				std::size_t v = set_;                                          \
				std::size_t c = 32;                                            \
				v &= -static_cast <std::make_signed<std::size_t>::type> (v);   \
                                                                               \
				if (v) c--;                                                    \
				if (v & 0x0000FFFF) c -= 16;                                   \
				if (v & 0x00FF00FF) c -= 8;                                    \
				if (v & 0x0F0F0F0F) c -= 4;                                    \
				if (v & 0x33333333) c -= 2;                                    \
				if (v & 0x55555555) c -= 1;                                    \
                                                                               \
				return static_cast <name> (c);                                 \
			}                                                                  \
                                                                               \
		private:                                                               \
			std::size_t set_;                                                  \
		};                                                                     \
                                                                               \
		const_iterator begin() const                                           \
		{                                                                      \
			return const_iterator(set_);                                       \
		}                                                                      \
                                                                               \
		const_iterator end() const                                             \
		{                                                                      \
			return const_iterator();                                           \
		}                                                                      \
                                                                               \
		name ## _set()                                                         \
			: set_{}                                                           \
		{                                                                      \
		}                                                                      \
                                                                               \
		name ## _set(std::initializer_list<name> list)                         \
			: set_(std::accumulate(list.begin(), list.end(), 0u,               \
				[](std::size_t set, name object) {                             \
					return set | (1 << static_cast <std::size_t> (object));    \
			  }))                                                              \
		{                                                                      \
		}                                                                      \
                                                                               \
		bool empty() const                                                     \
		{                                                                      \
			return set_ == 0;                                                  \
		}                                                                      \
                                                                               \
		void subtract(name ## _set other_set)                                  \
		{                                                                      \
			set_ ^= (set_ & other_set.set_);                                   \
		}                                                                      \
                                                                               \
		void add(name object)                                                  \
		{                                                                      \
			auto const index = static_cast <std::size_t> (object);             \
			set_ |= (1 << index);                                              \
		}                                                                      \
                                                                               \
		bool includes(name object) const                                       \
		{                                                                      \
			auto const index = static_cast <std::size_t> (object);             \
			return set_ & (1 << index);                                        \
		}                                                                      \
                                                                               \
		std::size_t size() const                                               \
		{                                                                      \
			return std::bitset<64>(set_).count();                              \
		}                                                                      \
                                                                               \
		bool operator==(name ## _set other) const                              \
		{                                                                      \
			return set_ == other.set_;                                         \
		}                                                                      \
                                                                               \
	private:                                                                   \
		std::size_t set_;                                                      \
	}

#endif
