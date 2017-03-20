#ifndef HEP_PS_PROTON_PDF_HPP
#define HEP_PS_PROTON_PDF_HPP

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

#include "hep/ps/parton.hpp"

#include <cstddef>
#include <memory>
#include <string>

namespace hep
{

/// Class that gives access to the proton pdfs.
template <typename T>
class proton_pdf
{
public:
	/// Constructor.
	proton_pdf(std::string const& name, std::size_t pdf_member);

	/// Move constructor.
	proton_pdf(proton_pdf<T>&& pdf);

	/// Destructor.
	~proton_pdf();

	/// Returns the strong coupling constant for the given `scale` compatible
	/// with the current parton distribution function.
	T alphas(T scale);

	/// Returns the value of the parton distribution function for the given \f$
	/// x \f$ and \f$ Q \f$, the last given by `scale`.
	parton_array<T> pdf(T x, T scale);

private:
	class impl;
	std::unique_ptr<impl> pimpl;
};

}

#endif
