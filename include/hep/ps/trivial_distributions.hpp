#ifndef HEP_PS_TRIVIAL_DISTRIBUTIONS_HPP
#define HEP_PS_TRIVIAL_DISTRIBUTIONS_HPP

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

#include "hep/mc/projector.hpp"

#include "hep/ps/cut_result.hpp"
#include "hep/ps/event_type.hpp"
#include "hep/ps/neg_pos_results.hpp"

#include <vector>

namespace hep
{

/// Trivial distribution generator, i.e. this functor does not generate any
/// distributions at all.
template <typename T>
class trivial_distributions
{
public:
	/// Constructor.
	trivial_distributions() = default;

	/// Does nothing.
	template <typename I>
	void operator()(
		std::vector<T> const&,
		cut_result_with_info<I> const&,
		neg_pos_results<T> const&,
		T,
		event_type,
		projector<T>&
	) {
	}

	/// Does nothing.
	template <typename I>
	void operator()(
		std::vector<T> const&,
		cut_result_with_info<I> const&,
		std::vector<neg_pos_results<T>> const&,
		T,
		event_type,
		projector<T>&
	) {
	}

	/// Does nothing.
	template <typename I>
	void operator()(
		std::vector<T> const&,
		cut_result_with_info<I> const&,
		std::vector<neg_pos_results<T>> const&,
		std::vector<neg_pos_results<T>> const&,
		T,
		event_type,
		projector<T>&
	) {
	}
};

}

#endif
