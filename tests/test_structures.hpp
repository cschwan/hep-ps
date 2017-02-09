#ifndef TEST_STRUCTURES_HPP
#define TEST_STRUCTURES_HPP

#include "hep/ps/cut_result.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_invariants.hpp"
#include "hep/ps/dipole_type.hpp"
#include "hep/ps/event_type.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/initial_state_array.hpp"
#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/particle_type.hpp"
#include "hep/ps/scales.hpp"

#include "catch.hpp"

#include <array>
#include <cstddef>
#include <vector>
#include <numeric>

namespace
{

template <typename T>
class test_luminosities
{
public:
	test_luminosities(hep::initial_state_set set)
		: set_(set)
	{
	}

	hep::initial_state_array<T> pdfs(T x1, T x2, T)
	{
		REQUIRE( x1 >= T() );
		REQUIRE( x1 < T(1.0) );
		REQUIRE( x2 >= T() );
		REQUIRE( x2 < T(1.0) );

		hep::initial_state_array<T> result;

		for (auto const process : hep::initial_state_list())
		{
			if (set_.includes(process))
			{
				result.set(process, T(1.0));
			}
		}

		return result;
	}

private:
	hep::initial_state_set set_;
};

template <typename T>
class test_matrix_elements
{
public:
	test_matrix_elements(
		hep::initial_state_set set,
		std::size_t final_states,
		std::size_t emitter,
		std::size_t unresolved,
		std::size_t spectator
	)
		: set_(set)
		, final_states_(final_states)
		, emitter_(emitter)
		, unresolved_(unresolved)
		, spectator_(spectator)
	{
	}

	T dipole(
		std::vector<T> const& dipole_phase_space,
		hep::initial_state process,
		hep::dipole const& dipole_info
	) const {
		REQUIRE( dipole_phase_space.size() == 4 * (final_states_ + 2) );
		REQUIRE( set_.includes(process) );
		REQUIRE( dipole_info.emitter() == emitter_ );
		REQUIRE( dipole_info.unresolved() == unresolved_ );
		REQUIRE( dipole_info.spectator() == spectator_ );

		return T(0.5);
	}

	std::vector<hep::dipole> dipole_ids(hep::initial_state state) const
	{
		REQUIRE( set_.includes(state) );

		return { hep::dipole(emitter_, unresolved_, spectator_,
			hep::particle_type::fermion, hep::particle_type::boson,
			hep::particle_type::fermion, hep::dipole_type::final_initial) };
	}

	hep::initial_state_array<T> borns(
		std::vector<T> const& phase_space,
		hep::initial_state_set set
	) const {
		REQUIRE( phase_space.size() == 4 * (final_states_ + 2) );
		REQUIRE( set == set_ );

		hep::initial_state_array<T> result;

		for (auto const process : hep::initial_state_list())
		{
			if (set.includes(process))
			{
				result.set(process, T(1.0));
			}
		}

		return result;
	}

	hep::initial_state_array<T> reals(
		std::vector<T> const& real_phase_space,
		hep::initial_state_set set
	) const {
		REQUIRE( real_phase_space.size() == 4 * (final_states_ + 3) );
		REQUIRE( set == set_ );

		hep::initial_state_array<T> result;

		for (auto const process : hep::initial_state_list())
		{
			if (set.includes(process))
			{
				result.set(process, T(1.0));
			}
		}

		return result;
	}

	std::vector<std::size_t> born_recombination_candidates() const
	{
		std::vector<std::size_t> indices(final_states_);
		std::iota(indices.begin(), indices.end(), 2);
		return indices;
	}

	std::vector<std::size_t> real_recombination_candidates() const
	{
		std::vector<std::size_t> indices(final_states_ + 1);
		std::iota(indices.begin(), indices.end(), 2);
		return indices;
	}

	void scale(T, test_luminosities<T>&)
	{
	}

private:
	hep::initial_state_set set_;
	std::size_t final_states_;
	std::size_t emitter_;
	std::size_t unresolved_;
	std::size_t spectator_;
};

template <typename T>
class test_subtraction
{
public:
	test_subtraction(
		std::size_t emitter,
		std::size_t unresolved,
		std::size_t spectator
	)
		: emitter_(emitter)
		, unresolved_(unresolved)
		, spectator_(spectator)
	{
	}

	hep::dipole_invariants<T> map_phase_space(
		std::vector<T> const& real_phase_space,
		std::vector<T>& dipole_phase_space,
		hep::dipole const& dipole_info
	) const {
		REQUIRE( (dipole_phase_space.size() + 4) == real_phase_space.size() );
		REQUIRE( dipole_info.emitter() == emitter_ );
		REQUIRE( dipole_info.unresolved() == unresolved_ );
		REQUIRE( dipole_info.spectator() == spectator_ );

		return hep::dipole_invariants<T>(1.0, 2.0, 4.0, 8.0);
	}

	T fermion_function(
		hep::dipole const& dipole_info,
		hep::dipole_invariants<T> const& invariants
	) const {
		REQUIRE( dipole_info.emitter() == emitter_ );
		REQUIRE( dipole_info.unresolved() == unresolved_ );
		REQUIRE( dipole_info.spectator() == spectator_ );

		REQUIRE( invariants.one == T(1.0) );
		REQUIRE( invariants.two == T(2.0) );
		REQUIRE( invariants.sij == T(4.0) );
		REQUIRE( invariants.adipole == T(8.0) );

		return T(1.0);
	}

private:
	std::size_t emitter_;
	std::size_t unresolved_;
	std::size_t spectator_;
};

template <typename T>
struct test_cuts
{
public:
	test_cuts(std::size_t final_states, std::size_t inclusive)
		: final_states_(final_states)
		, inclusive_(inclusive)
	{
	}

	hep::cut_result cut(
		std::vector<T> const& phase_space,
		T,
		hep::event_type type
	) const {
		if (type == hep::event_type::inclusive_n_plus_1)
		{
			REQUIRE( phase_space.size() == 4 * (final_states_ + 3) );
		}
		else if (type == hep::event_type::born_like_n)
		{
			REQUIRE( phase_space.size() == 4 * (final_states_ + 2) );
		}
		else
		{
			FAIL( "type is neither born nor inclusive" );
		}

		if (type == hep::event_type::inclusive_n_plus_1 && !inclusive_)
		{
			return { true, true };
		}

		return { false, false };
	}

private:
	std::size_t final_states_;
	std::size_t inclusive_;
};

template <typename T>
class test_recombiner
{
public:
	test_recombiner(std::size_t final_states)
		: final_states_(final_states)
	{
	}

	std::size_t recombine(
		std::vector<T> const& phase_space,
		std::vector<T>& recombined_phase_space,
		std::vector<std::size_t> const& recombination_candidates,
		std::size_t
	) const {
		REQUIRE( phase_space.size() == recombined_phase_space.size() );
		REQUIRE( recombination_candidates.size() ==
			(phase_space.size() / 4 - 2) );

		// TODO: improve the test

		std::size_t last_index = 1;
		for (auto const index : recombination_candidates)
		{
			REQUIRE( index > last_index );
			last_index = index;
		}

		return 0;
	}

private:
	std::size_t final_states_;
};

template <typename T>
class test_scale_setter
{
public:
	static constexpr T scale = T(1.0);

	hep::scales<T> operator()(std::vector<T> const&)
	{
		return { scale , scale };
	}
};

}

#endif
