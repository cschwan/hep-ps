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
#include <cstddef>
#include <initializer_list>

// TODO: make the `enum` an `enum class`

/// Defines an enumeration `name`, a `constexpr` function `name_list()`, and a
/// template class `name_array<T>`. The arguments for this macro after `name`
/// are the possible values of the enumeration. The function returns a
/// `std::initializer_list` with the possible values of the enumeration. The
/// class `name_array<T>` is a wrapper over `std::array<T, ...>` that supports
/// access to its values by keys that are values of the enumeration.
#define HEP_ENUM(name, ...)                                                    \
	enum name : std::size_t                                                    \
	{                                                                          \
		__VA_ARGS__                                                            \
	};                                                                         \
	constexpr std::initializer_list<name> name ## _list()                      \
	{                                                                          \
		return { __VA_ARGS__ };                                                \
	}                                                                          \
	template <typename T>                                                      \
	class name ## _array                                                       \
	{                                                                          \
	public:                                                                    \
		name ## _array()                                                       \
			: array_{{}}                                                       \
		{                                                                      \
		}                                                                      \
	                                                                           \
		T& operator[](name index)                                              \
		{                                                                      \
			return array_[static_cast <std::size_t> (index)];                  \
		}                                                                      \
	                                                                           \
		T operator[](name index) const                                         \
		{                                                                      \
			return array_[static_cast <std::size_t> (index)];                  \
		}                                                                      \
	private:                                                                   \
		std::array<T, name ## _list().size()> array_;                          \
	}

#endif
