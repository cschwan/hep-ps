#ifndef HEP_PS_ALPHAS_CALC_HPP
#define HEP_PS_ALPHAS_CALC_HPP

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

namespace hep
{

/// Calculates the strong coupling at different scales.
template <typename T>
class alphas_calc
{
public:
	/// Constructor. Do not call, use of member functions of PDF classes
	/// instead.
	alphas_calc(void* pdf);

	alphas_calc(alphas_calc<T> const&) = delete;

	alphas_calc(alphas_calc<T>&&) = delete;

	alphas_calc& operator=(alphas_calc<T> const&) = delete;

	alphas_calc& operator=(alphas_calc<T>&&) = delete;

	/// Returns the strong coupling for the given `scale`.
	T alphas(T scale);

private:
	void* pdf_;
};

}

#endif
