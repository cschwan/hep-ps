#ifndef HEP_PS_LUMINOSITY_INFO_HPP
#define HEP_PS_LUMINOSITY_INFO_HPP

/*
 * hep-ps - A C++ Library for Perturbative Calculations in High Energy Physics
 * Copyright (C) 2016  Christopher Schwan
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

namespace hep
{

template <typename T>
class luminosity_info
{
public:
	luminosity_info(T x1, T x2, T energy_squared, T rapidity_shift)
		: x1_(x1)
		, x2_(x2)
		, energy_squared_(energy_squared)
		, rapidity_shift_(rapidity_shift)
	{
	}

	T x1() const
	{
		return x1_;
	}

	T x2() const
	{
		return x2_;
	}

	T energy_squared() const
	{
		return energy_squared_;
	}

	T rapidity_shift() const
	{
		return rapidity_shift_;
	}

private:
	T x1_;
	T x2_;
	T energy_squared_;
	T rapidity_shift_;
};

}

#endif
