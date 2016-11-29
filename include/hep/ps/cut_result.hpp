#ifndef HEP_PS_CUT_RESULT_HPP
#define HEP_PS_CUT_RESULT_HPP

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

class cut_result
{
public:
	cut_result(bool neg_cutted, bool pos_cutted)
		: neg_cutted_(neg_cutted)
		, pos_cutted_(pos_cutted)
	{
	}

	bool neg_cutted() const
	{
		return neg_cutted_;
	}

	bool pos_cutted() const
	{
		return pos_cutted_;
	}

private:
	bool neg_cutted_;
	bool pos_cutted_;
};

}

#endif
