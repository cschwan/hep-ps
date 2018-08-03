#ifndef HEP_PS_INSERTION_TERM_HPP
#define HEP_PS_INSERTION_TERM_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017-2018  Christopher Schwan
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

#include "hep/ps/correction_type.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/particle_type.hpp"

#include <cassert>
#include <cstddef>

namespace hep
{

/// Classifies a finite insertion term according to the emitter and spectator
/// particle, similarly to \ref dipole.
class insertion_term
{
public:
	/// Constructor for insertion terms not of the type `born`.
	insertion_term(
		std::size_t emitter,
		particle_type emitter_type,
		std::size_t spectator,
		correction_type corr_type
	)
		: emitter_(emitter)
		, emitter_type_(emitter_type)
		, spectator_(spectator)
		, corr_type_(corr_type)
	{
		int type = (emitter < 2) | ((spectator < 2) << 1);

		switch (type)
		{
		case 0:
			type_ = insertion_term_type::final_final;
			break;

		case 1:
			type_ = insertion_term_type::initial_final;
			break;

		case 2:
			type_ = insertion_term_type::final_initial;
			break;

		case 3:
			type_ = insertion_term_type::initial_initial;
			break;
		}
	}

	/// Constructor for `born` insertion terms.
	insertion_term(std::size_t initial, correction_type corr_type)
		: type_(insertion_term_type::born)
		, emitter_{initial}
		, emitter_type_{particle_type::boson}
		, spectator_{0}
		, corr_type_(corr_type)
	{
	}

	/// Returns the index of the emitter particle.
	std::size_t emitter() const
	{
		assert( type_ != insertion_term_type::born );

		return emitter_;
	}

	/// Returns the index of the spectator particle.
	std::size_t spectator() const
	{
		assert( type_ != insertion_term_type::born );

		return spectator_;
	}

	/// Returns the type of the emitter particle.
	particle_type emitter_type() const
	{
		assert( type_ != insertion_term_type::born );

		return emitter_type_;
	}

	/// Returns the type of this insertion term.
	insertion_term_type type() const
	{
		return type_;
	}

	/// Returns the index of the initial particle.
	std::size_t initial_particle() const
	{
		switch (type_)
		{
		case insertion_term_type::born:
			return emitter_;

		case insertion_term_type::initial_final:
		case insertion_term_type::initial_initial:
			return emitter();

		case insertion_term_type::final_initial:
			return spectator();

		case insertion_term_type::final_final:
			assert( false );

		default:
			assert( false );
		}
	}

	/// Returns whether this insertion term is a QCD or an EW one.
	correction_type corr_type() const
	{
		return corr_type_;
	}

private:
	insertion_term_type type_;
	std::size_t emitter_;
	particle_type emitter_type_;
	std::size_t spectator_;
	correction_type corr_type_;
};

/// Comparison operator for two instances of \ref insertion_term `a` and `b`.
bool operator<(insertion_term const& a, insertion_term const& b);

/// Comparison operator for two instances of \ref insertion_term `a` and `b`.
bool operator==(insertion_term const& a, insertion_term const& b);

}

namespace std
{

template <>
struct hash<hep::insertion_term>
{
	size_t operator()(hep::insertion_term const& term) const;
};

}

#endif
