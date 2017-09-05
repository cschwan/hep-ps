#ifndef TEST_STRUCTURES_HPP
#define TEST_STRUCTURES_HPP

#include "hep/ps/abc_terms.hpp"
#include "hep/ps/cut_result.hpp"
#include "hep/ps/dipole.hpp"
#include "hep/ps/dipole_with_set.hpp"
#include "hep/ps/dipole_invariants.hpp"
#include "hep/ps/dipole_type.hpp"
#include "hep/ps/event_type.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/insertion_term.hpp"
#include "hep/ps/insertion_term_type.hpp"
#include "hep/ps/particle_type.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/scales.hpp"

#include "catch.hpp"

#include <array>
#include <cstddef>
#include <vector>
#include <numeric>

namespace
{

template <typename T>
class test_pdf
{
public:
	hep::parton_array<T> pdf(T x, T scale)
	{
		CHECK( x >= T() );
		CHECK( x < T(1.0) );
		CHECK( scale > T() );

		hep::parton_array<T> result;

		result[hep::parton::anti_up]      = T(1.0);
		result[hep::parton::anti_down]    = T(1.0);
		result[hep::parton::anti_charm]   = T(1.0);
		result[hep::parton::anti_strange] = T(1.0);
		result[hep::parton::gluon]        = T(1.0);
		result[hep::parton::up]           = T(1.0);
		result[hep::parton::down]         = T(1.0);
		result[hep::parton::charm]        = T(1.0);
		result[hep::parton::strange]      = T(1.0);

		return result;
	}

	int alphas()
	{
		return 0;
	}
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

	hep::initial_state_array<T> dipole_me(
		hep::dipole const& dipole_info,
		std::vector<T> const& dipole_phase_space,
		hep::initial_state_set set
	) const {
		CHECK( dipole_phase_space.size() == 4 * (final_states_ + 2) );
		CHECK( set_ == set );
		CHECK( dipole_info.emitter() == emitter_ );
		CHECK( dipole_info.unresolved() == unresolved_ );
		CHECK( dipole_info.spectator() == spectator_ );

		hep::initial_state_array<T> result;

		for (auto const state : set_)
		{
			result[state] = T(0.5);
		}

		return result;
	}

	std::array<hep::dipole_with_set, 1> dipoles() const
	{
		return { hep::dipole_with_set(emitter_, unresolved_, spectator_,
			hep::particle_type::fermion, hep::particle_type::boson,
			hep::particle_type::fermion, hep::dipole_type::final_initial, set_)
		};
	}

	hep::initial_state_array<T> borns(
		std::vector<T> const& phase_space,
		hep::initial_state_set set
	) const {
		CHECK( phase_space.size() == 4 * (final_states_ + 2) );
		CHECK( set == set_ );

		hep::initial_state_array<T> result;

		for (auto const process : hep::initial_state_list())
		{
			if (set.includes(process))
			{
				result[process] = T(1.0);
			}
		}

		return result;
	}

	hep::initial_state_array<T> reals(
		std::vector<T> const& real_phase_space,
		hep::initial_state_set set
	) const {
		CHECK( real_phase_space.size() == 4 * (final_states_ + 3) );
		CHECK( set == set_ );

		hep::initial_state_array<T> result;

		for (auto const process : hep::initial_state_list())
		{
			if (set.includes(process))
			{
				result[process] = T(1.0);
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

	std::array<hep::insertion_term, 3> insertion_terms() const
	{
		return {
			hep::insertion_term(
				hep::insertion_term_type::born
			),
			hep::insertion_term(
				hep::insertion_term_type::final_initial,
				emitter_,
				hep::particle_type::fermion,
				spectator_
			),
			hep::insertion_term(
				hep::insertion_term_type::initial_final,
				emitter_,
				hep::particle_type::fermion,
				spectator_
			)
		};
	}

	std::array<hep::initial_state_array<T>, 3> correlated_me(
		std::vector<T> const&,
		hep::initial_state_set
	) {
		return {
			hep::initial_state_array<T>(),
			hep::initial_state_array<T>(),
			hep::initial_state_array<T>()
		};
	}

	void scale(T, int)
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
		CHECK( (dipole_phase_space.size() + 4) == real_phase_space.size() );
		CHECK( dipole_info.emitter() == emitter_ );
		CHECK( dipole_info.unresolved() == unresolved_ );
		CHECK( dipole_info.spectator() == spectator_ );

		return hep::dipole_invariants<T>(1.0, 2.0, 4.0, 8.0);
	}

	T fermion_function(
		hep::dipole const& dipole_info,
		hep::dipole_invariants<T> const& invariants
	) const {
		CHECK( dipole_info.emitter() == emitter_ );
		CHECK( dipole_info.unresolved() == unresolved_ );
		CHECK( dipole_info.spectator() == spectator_ );

		CHECK( invariants.one == T(1.0) );
		CHECK( invariants.two == T(2.0) );
		CHECK( invariants.sij == T(4.0) );
		CHECK( invariants.alpha == T(8.0) );

		return T(1.0);
	}

	hep::abc_terms<T> insertion_terms(
		hep::insertion_term const&,
		hep::scales<T> const&,
		std::vector<T> const&,
		T,
		T
	) const {
		return hep::abc_terms<T>();
	}

	T insertion_terms2(
		hep::insertion_term const&,
		hep::scales<T> const&,
		std::vector<T> const&
	) const {
		return T();
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
			CHECK( phase_space.size() == 4 * (final_states_ + 3) );
		}
		else if (type == hep::event_type::born_like_n)
		{
			// FIXME: this check fails for finite insertion terms
//			CHECK( phase_space.size() == 4 * (final_states_ + 2) );
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
		CHECK( phase_space.size() == recombined_phase_space.size() );
		CHECK( recombination_candidates.size() ==
			(phase_space.size() / 4 - 2) );

		// TODO: improve the test

		std::size_t last_index = 1;
		for (auto const index : recombination_candidates)
		{
			CHECK( index > last_index );
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

template <typename T>
inline std::vector<hep::distribution_parameters<T>>
test_distribution_parameters()
{
	return std::vector<hep::distribution_parameters<T>>();
}

}

#endif
