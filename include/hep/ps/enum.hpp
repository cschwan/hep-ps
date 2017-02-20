#ifndef HEP_PS_ENUM_HPP
#define HEP_PS_ENUM_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
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

/// Defines an enumeration `name`. The arguments after `name` are the possible
/// values of the enumeration.
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
	using name ## _array_ = std::array<T, name ## _list().size()>

#endif
