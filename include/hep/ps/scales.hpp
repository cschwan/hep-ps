#ifndef HEP_PS_SCALES_HPP
#define HEP_PS_SCALES_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
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

namespace hep
{

/// Class that captures the factorization and renormalization scales.
template <typename T>
class scales
{
public:
	/// Constructs a new object with \f$ \mu_\text{F} \f$ set to the value of
	/// `factorization` and \f$ \mu_\text{R} \f$ set to the value of
	/// `renormalization`.
	scales(T factorization, T renormalization)
		: factorization_(factorization)
		, renormalization_(renormalization)
	{
	}

	/// Returns the factorization scale \f$ \mu_\text{F} \f$.
	T factorization() const
	{
		return factorization_;
	}

	/// Returns the renormalization scale \f$ \mu_\text{R} \f$.
	T renormalization() const
	{
		return renormalization_;
	}

private:
	T factorization_;
	T renormalization_;
};

}

#endif
