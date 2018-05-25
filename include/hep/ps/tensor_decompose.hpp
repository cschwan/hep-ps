#ifndef HEP_PS_TENSOR_DECOMPOSE_HPP
#define HEP_PS_TENSOR_DECOMPOSE_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2018  Christopher Schwan
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

namespace hep
{

/// Decomposes a tensor \f$ T^{\mu \nu} \f$ into a sum of product of vectors:
/// \f[
///   T^{\mu \nu} = - g^{\mu \nu} + \frac{p^\mu p^\nu}{D p^2} =
///   \sum_{\lambda \in \{ 1, 2, 3, 4 \} } f_\lambda \varepsilon_\lambda^{\mu}
///   \varepsilon_\lambda^{\nu}
/// \f]
/// where the variable \f$ D \f$ is given by `denominator`, the momentum \f$ p
/// \f$ by `momentum`. The resulting vectors \f$ \varepsilon_\lambda \f$ are
/// written into `results` and the variables \f$ f_\lambda \f$ are written into
/// `factors`.
template <typename T>
void tensor_decompose(
	std::array<T, 4> const& momentum,
	T denominator,
	std::array<std::array<T, 4>, 4>& results,
	std::array<T, 4>& factors
);

}

#endif
