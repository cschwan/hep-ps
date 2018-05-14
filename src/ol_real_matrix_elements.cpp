#include "hep/ps/ol_interface.hpp"
#include "hep/ps/ol_real_matrix_elements.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>
#include <utility>

namespace hep
{

template <typename T>
ol_real_matrix_elements<T>::ol_real_matrix_elements(
	std::vector<std::string> const& real_processes,
	std::size_t alphas_power,
	correction_type type
)
	: alphas_power_(alphas_power)
	, type_(type)
{
	auto& ol = ol_interface::instance();

	std::size_t const dipole_qcd_order = (type == correction_type::qcd) ?
		(alphas_power - 1) : alphas_power;

	std::multimap<dipole, initial_state> dipoles_with_state;

	for (auto const& process : real_processes)
	{
		auto const& pdg_ids = ol_process_string_to_pdg_ids(process);

		std::vector<final_state> final_states(pdg_ids.size() - 2);
		std::transform(pdg_ids.begin() + 2, pdg_ids.end(), final_states.begin(),
			pdg_id_to_final_state);

		auto const state = partons_to_initial_state(
			pdg_id_to_parton(pdg_ids.at(0)), pdg_id_to_parton(pdg_ids.at(1)));

		if (final_states_real_.empty())
		{
			final_states_real_ = final_states;
		}
		else
		{
			if (!std::equal(final_states.begin(), final_states.end(),
				final_states_real_.begin()))
			{
				throw std::invalid_argument("processes are not compatible");
			}
		}

		ol.setparameter_int("order_qcd", alphas_power);
		ids_reals_.emplace(state, ol.register_process(process.c_str(), 1));

		std::vector<std::size_t> unresolved;
		std::vector<std::size_t> indices;

		if (type == correction_type::ew)
		{
			// TODO: NYI
			assert( false );

			// TODO: initialize `charge_table`
		}
		else if (type == correction_type::qcd)
		{
			// TODO: the following construction only works if there are no
			// initial state gluons and if there are no gluon -> fermion
			// anti-fermion splittings

			std::vector<final_state> final_states;

			for (std::size_t i = 0; i != pdg_ids.size(); ++i)
			{
				if (pdg_ids.at(i) == 21)
				{
					unresolved.push_back(i);
				}
				else if (pdg_id_particle_has_color(pdg_ids.at(i)))
				{
					indices.push_back(i);
				}
			}

			for (auto const un : unresolved)
			{
				final_states = final_states_real_;
				final_states.erase(final_states.begin() + un - 2);

				if (final_states_.empty())
				{
					final_states_ = final_states;
				}
				else
				{
					assert( std::equal(final_states.begin(), final_states.end(),
						final_states_.begin()) );
				}

				std::string dipole_process;
				dipole_process.append(std::to_string(pdg_ids.at(0)));
				dipole_process.append(" ");
				dipole_process.append(std::to_string(pdg_ids.at(1)));
				dipole_process.append(" -> ");

				for (std::size_t i = 2; i != pdg_ids.size(); ++i)
				{
					if (i != un)
					{
						dipole_process.append(std::to_string(pdg_ids.at(i)));
						dipole_process.append(" ");
					}
				}

				ol.setparameter_int("order_qcd", dipole_qcd_order);

				ids_dipoles_.emplace(state,
					ol.register_process(dipole_process.c_str(), 1));

				auto const un_t = pdg_id_to_particle_type(pdg_ids.at(un));

				for (std::size_t i = 0; i < (indices.size() - 1); ++i)
				{
					std::size_t const em = indices.at(i);
					auto const em_t = pdg_id_to_particle_type(pdg_ids.at(em));

					for (std::size_t k = i + 1; k != indices.size(); ++k)
					{
						std::size_t const sp = indices.at(k);
						auto const sp_t = pdg_id_to_particle_type(
							pdg_ids.at(sp));

						dipoles_with_state.insert(std::make_pair(dipole(em, un,
							sp, em_t, un_t, sp_t), state));
						dipoles_with_state.insert(std::make_pair(dipole(sp, un,
							em, sp_t, un_t, em_t), state));
					}
				}
			}
		}
		else
		{
			assert( false );
		}
	}

	std::vector<dipole> dipoles;
	for (auto const& dip : dipoles_with_state)
	{
		auto const& dipole = dip.first;

		if (std::find(dipoles.begin(), dipoles.end(), dipole) == dipoles.end())
		{
			dipoles.push_back(dipole);
		}
	}

	for (auto const& dipole : dipoles)
	{
		initial_state_set set;

		auto const range = dipoles_with_state.equal_range(dipole);

		for (auto i = range.first; i != range.second; ++i)
		{
			set.add(i->second);
		}

		dipoles_.emplace_back(
			dipole.emitter(),
			dipole.unresolved(),
			dipole.spectator(),
			dipole.emitter_type(),
			dipole.unresolved_type(),
			dipole.spectator_type(),
			set
		);
	}

	// TODO: merge entries of `dipoles_` which have the same dipole but a
	// different set

	final_states_real_.shrink_to_fit();
	final_states_.shrink_to_fit();
	std::size_t const n = final_states_.size() + 2;
	ol_m2cc_.resize(n * (n - 1) / 2);
	ol_phase_space_.resize(5 * (n + 1));
}

template <typename T>
void ol_real_matrix_elements<T>::alphas(T alphas)
{
	auto& ol = ol_interface::instance();
	ol.setparameter_double("alphas", static_cast <double>(alphas));
}

template <typename T>
std::size_t ol_real_matrix_elements<T>::alphas_power() const
{
	return alphas_power_;
}

template <typename T>
void ol_real_matrix_elements<T>::dipole_me(
	dipole const& dipole,
	std::vector<T> const& phase_space,
	initial_state_set set,
	std::vector<scales<T>> const& scales,
	std::vector<initial_state_map<T>>& results
) {
	auto& ol = hep::ol_interface::instance();

	// TODO: for the time being we assume that the matrix elements do only
	// indirectly depend on the renormalization scales (through alphas)

	std::size_t const n = phase_space.size() / 4;

	for (std::size_t i = 0; i != n; ++i)
	{
		ol_phase_space_.at(5 * i + 0) =
			static_cast <double> (phase_space.at(4 * i + 0));
		ol_phase_space_.at(5 * i + 1) =
			static_cast <double> (phase_space.at(4 * i + 1));
		ol_phase_space_.at(5 * i + 2) =
			static_cast <double> (phase_space.at(4 * i + 2));
		ol_phase_space_.at(5 * i + 3) =
			static_cast <double> (phase_space.at(4 * i + 3));
		ol_phase_space_.at(5 * i + 4) = 0.0;
	}

	if (type_ == correction_type::qcd)
	{
		double m2tree;
		double m2ew;
		double alphas;
		ol.getparameter_double("alphas", &alphas);

		for (auto const state : set)
		{
			auto const range = ids_dipoles_.equal_range(state);

			for (auto i = range.first; i != range.second; ++i)
			{
				ol.evaluate_cc(i->second, ol_phase_space_.data(), &m2tree,
					ol_m2cc_.data(), &m2ew);

				std::size_t const k = std::min(dipole.emitter(),
					dipole.spectator());
				std::size_t const l = std::max(dipole.emitter(),
					dipole.spectator());
				std::size_t const index = k + l * (l - 1) / 2;

				results.front().emplace_back(state,
					T(alphas) * T(ol_m2cc_.at(index)));
			}
		}
	}
	else if (type_ == correction_type::ew)
	{
		double m2tree;
		double alpha;
		ol.getparameter_double("alpha", &alpha);

		for (auto const state : set)
		{
			auto const range = ids_dipoles_.equal_range(state);

			for (auto i = range.first; i != range.second; ++i)
			{
				ol.evaluate_tree(i->second, ol_phase_space_.data(), &m2tree);

				T const charge_em = charge_table_.at(i->second).at(
					dipole.emitter());
				T const charge_sp = charge_table_.at(i->second).at(
					dipole.spectator());

				results.front().emplace_back(state,
					charge_em * charge_sp * T(alpha) * T(m2tree));
			}
		}
	}
	else
	{
		assert( false );
	}

	for (std::size_t i = 1; i != scales.size(); ++i)
	{
		results.at(i) = results.front();
	}
}

template <typename T>
std::vector<dipole_with_set> const&
ol_real_matrix_elements<T>::dipoles() const
{
	return dipoles_;
}

template <typename T>
std::vector<final_state> const&
ol_real_matrix_elements<T>::final_states() const
{
	return final_states_;
}

template <typename T>
std::vector<final_state> const&
ol_real_matrix_elements<T>::final_states_real() const
{
	return final_states_real_;
}

template <typename T>
void ol_real_matrix_elements<T>::reals(
	std::vector<T> const& phase_space,
	initial_state_set set,
	std::vector<scales<T>> const& scales,
	std::vector<initial_state_map<T>>& results
) {
	auto& ol = hep::ol_interface::instance();

	// TODO: for the time being we assume that the matrix elements do only
	// indirectly depend on the renormalization scales (through alphas)

	std::size_t const n = phase_space.size() / 4;

	for (std::size_t i = 0; i != n; ++i)
	{
		ol_phase_space_.at(5 * i + 0) =
			static_cast <double> (phase_space.at(4 * i + 0));
		ol_phase_space_.at(5 * i + 1) =
			static_cast <double> (phase_space.at(4 * i + 1));
		ol_phase_space_.at(5 * i + 2) =
			static_cast <double> (phase_space.at(4 * i + 2));
		ol_phase_space_.at(5 * i + 3) =
			static_cast <double> (phase_space.at(4 * i + 3));
		ol_phase_space_.at(5 * i + 4) = 0.0;
	}

	double m2tree;

	for (auto const state : set)
	{
		auto const range = ids_reals_.equal_range(state);
		T value = T();

		for (auto i = range.first; i != range.second; ++i)
		{
			ol.evaluate_tree(i->second, ol_phase_space_.data(), &m2tree);
			value += T(m2tree);
		}

		results.front().emplace_back(state, value);
	}

	for (std::size_t i = 1; i != scales.size(); ++i)
	{
		results.at(i) = results.front();
	}
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class ol_real_matrix_elements<double>;

}
