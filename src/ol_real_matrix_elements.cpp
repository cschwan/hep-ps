#include "hep/ps/correction_type.hpp"
#include "hep/ps/ol_interface.hpp"
#include "hep/ps/ol_real_matrix_elements.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>
#include <utility>

namespace
{

std::vector<int> dipole_me(
	std::vector<int> const& pdg_ids,
	std::size_t i,
	std::size_t j,
	std::size_t k,
	hep::correction_type type,
	hep::coupling_order real_me_order
) {
	std::vector<int> result;

	int const id_i = pdg_ids.at(i);
	int const id_j = pdg_ids.at(j);
	int const id_k = pdg_ids.at(k);
	int const sign = ((i < 2) == (j < 2)) ? 1 : -1;

	if (type == hep::correction_type::ew)
	{
		if (hep::pdg_id_has_charge(id_i) && hep::pdg_id_has_charge(id_k))
		{
			if (id_j == hep::pdg_id_of_photon())
			{
				result = pdg_ids;
				result.erase(result.begin() + j);
			}
			else if (hep::pdg_id_has_charge(id_j) &&
				((id_i + sign * id_j) == 0))
			{
				result = pdg_ids;
				result.erase(result.begin() + j);
				result.at(i) = hep::pdg_id_of_photon();
			}
		}
	}
	else if (type == hep::correction_type::qcd)
	{
		if (hep::pdg_id_has_color(id_i) && hep::pdg_id_has_color(id_j) &&
			hep::pdg_id_has_color(id_k))
		{
			if (id_j == hep::pdg_id_of_gluon())
			{
				result = pdg_ids;
				result.erase(result.begin() + j);
			}
			else if ((id_i + sign * id_j) == 0)
			{
				result = pdg_ids;
				result.erase(result.begin() + j);
				result.at(i) = hep::pdg_id_of_gluon();
			}
		}
	}
	else
	{
		assert( false );
	}

	std::size_t const gluons = std::count_if(result.begin(), result.end(),
		hep::pdg_id_is_gluon);
	std::size_t const quarks = std::count_if(result.begin(), result.end(),
		hep::pdg_id_is_quark);

	std::size_t n_qcd_min;
	std::size_t n_qcd_max;

	if (gluons >= (result.size() - 2))
	{
		// pure gluon case and one quark line
		n_qcd_min = (gluons == result.size()) ? (gluons - 2) : gluons;
		n_qcd_max = n_qcd_min;
	}
	else
	{
		// starting with two quark lines it is possible to have photons or
		// gluons in between the quark lines
		n_qcd_min = gluons;
		n_qcd_max = n_qcd_min + quarks / 2;
	}

	std::size_t order_qcd = (type == hep::correction_type::qcd) ?
		(real_me_order.alphas_power() - 1) : real_me_order.alphas_power();

	// check if the matrix element with given coupling order exists
	if ((order_qcd < n_qcd_min) || (order_qcd > n_qcd_max))
	{
		result.clear();
	}

	return result;
}

}

namespace hep
{

template <typename T>
ol_real_matrix_elements<T>::ol_real_matrix_elements(
	std::vector<std::string> const& real_processes,
	std::vector<final_state> const& dipole_final_states,
	coupling_order order
)
	: alphas_power_(order.alphas_power())
	, final_states_(dipole_final_states)
{
	auto& ol = ol_interface::instance();

	for (auto const& process : real_processes)
	{
		auto const& pdg_ids = ol_process_string_to_pdg_ids(process);
		auto const& states = pdg_ids_to_states(pdg_ids);

		if (final_states_real_.empty())
		{
			final_states_real_ = states.second;
		}
		else
		{
			if (!std::equal(states.second.begin(), states.second.end(),
				final_states_real_.begin()))
			{
				throw std::invalid_argument("processes are not compatible");
			}
		}

		ol.setparameter_int("order_qcd", order.alphas_power());
		ids_reals_.emplace(states.first,
			ol.register_process(process.c_str(), 1));

		// construct all possible QCD and EW dipoles
		for (std::size_t i = 0; i != pdg_ids.size(); ++i)
		{
		for (std::size_t j = 2; j != pdg_ids.size(); ++j)
		{
		for (std::size_t k = 0; k != pdg_ids.size(); ++k)
		{
			if ((i == j) || (i == k) || (j == k))
			{
				continue;
			}

			for (auto const type : correction_type_list())
			{
				auto const& ids = ::dipole_me(pdg_ids, i, j, k, type, order);

				if (ids.empty())
				{
					// there is no matrix element for this dipole
					continue;
				}

				auto const& states = pdg_ids_to_states(ids);

				// if the signature does not match, throw the dipole away
				if (!std::equal(states.second.begin(), states.second.end(),
					final_states_.begin()))
				{
					continue;
				}

				ol.setparameter_int("order_qcd", (type == correction_type::qcd)
					? (order.alphas_power() - 1) : order.alphas_power());
				auto const process = pdg_ids_to_ol_process_string(ids);
				int const dipole_id = ol.register_process(process.c_str(), 1);

				if (type == correction_type::ew)
				{
					std::vector<T> charges;
					charges.reserve(ids.size());

					for (int const id : ids)
					{
						charges.push_back(T(pdg_id_to_charge_times_three(id)) /
							T(3.0));
					}

					charge_table_.emplace(dipole_id, std::move(charges));
				}

				ids_dipoles_.emplace(states.first, dipole_id);

				auto const type_i = pdg_id_to_particle_type(pdg_ids.at(i));
				auto const type_j = pdg_id_to_particle_type(pdg_ids.at(j));
				auto const type_k = pdg_id_to_particle_type(pdg_ids.at(k));

				dipoles_.emplace_back(i, j, k, type_i, type_j, type_k, type);
			}
		}
		}
		}
	}

	std::sort(dipoles_.begin(), dipoles_.end());
	auto new_end = std::unique(dipoles_.begin(), dipoles_.end());
	dipoles_.erase(new_end, dipoles_.end());
	dipoles_.shrink_to_fit();

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

	if (dipole.corr_type() == correction_type::qcd)
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
	else if (dipole.corr_type() == correction_type::ew)
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
std::vector<dipole> const& ol_real_matrix_elements<T>::dipoles() const
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
