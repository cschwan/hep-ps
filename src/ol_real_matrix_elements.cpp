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

	// TODO: technically we should also allow dipoles with uncharged particles
	// as spectators if they are proportional to EW spin-correlated dipoles

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
	coupling_order order,
	photon_dipole_selector const& selector
)
	: alphas_power_(order.alphas_power())
	, final_states_(dipole_final_states)
{
	auto& ol = ol_interface::instance();

	for (auto const& process : real_processes)
	{
		auto const& ids = ol_process_string_to_pdg_ids(process);
		auto const& states = pdg_ids_to_states(ids);

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

		// construct all possible EW and QCD dipoles
		for (auto const type : correction_type_list())
		{
		for (std::size_t i = 0; i != ids.size(); ++i)
		{
		for (std::size_t j = 2; j != ids.size(); ++j)
		{
		for (std::size_t k = 0; k != ids.size(); ++k)
		{
			if ((i == j) || (i == k) || (j == k))
			{
				continue;
			}

			auto const& dipole_ids = ::dipole_me(ids, i, j, k, type, order);

			if (dipole_ids.empty())
			{
				// there is no matrix element for this dipole
				continue;
			}

			auto const& dipole_states = pdg_ids_to_states(dipole_ids);

			auto const mismatch = std::mismatch(final_states_.begin(),
				final_states_.end(), dipole_states.second.begin());

			if (mismatch.first != final_states_.end())
			{
				// check if the mismatch is caused by a quark-antiquark-pair
				// that got recombined into a photon; if this is the case allow
				// the dipole
				if ((*mismatch.first != final_state::quark_gluon) ||
					(*mismatch.second != final_state::photon) ||
					!std::equal(mismatch.first + 1, final_states_.end(),
						mismatch.second + 1))
				{
					continue;
				}
			}

			bool photon_dipole_selected = false;

			if (type == correction_type::ew)
			{
				int const id_i = ids.at(i);
				int const id_j = ids.at(j);
				int const sign = ((i < 2) == (j < 2)) ? 1 : -1;

				// check if this is a photon dipole
				if (id_i + sign * id_j == 0)
				{
					if (selector(dipole_ids, i, j, k))
					{
						photon_dipole_selected = true;
					}
					else
					{
						continue;
					}
				}
			}

			ol.setparameter_int("order_qcd", (type == correction_type::qcd)
				? (order.alphas_power() - 1) : order.alphas_power());
			auto const process = pdg_ids_to_ol_process_string(dipole_ids);
			int const dipole_id = ol.register_process(process.c_str(), 1);

			// add a charge table if neccessary
			if ((type == correction_type::ew) &&
				(charge_table_.find(dipole_id) == charge_table_.end()))
			{
				std::vector<T> charges;
				charges.reserve(ids.size());

				// we need the charges of the real matrix element
				for (int const id : ids)
				{
					int const charge = pdg_id_to_charge_times_three(id);
					charges.push_back(T(charge) / T(3.0));
				}

				// sign of the crossing
				charges.at(0) *= T(-1.0);
				charges.at(1) *= T(-1.0);

				charge_table_.emplace(dipole_id, std::move(charges));
			}

			auto const type_i = pdg_id_to_particle_type(ids.at(i));
			auto const type_j = pdg_id_to_particle_type(ids.at(j));
			auto const type_k = pdg_id_to_particle_type(ids.at(k));

			auto const dip = dipole(i, j, k, type_i, type_j, type_k, type);

			// add a dipole if it doesn't exist yet
			if (std::find(dipoles_.begin(), dipoles_.end(), dip) ==
				dipoles_.end())
			{
				dipoles_.push_back(dip);
			}

			auto range = mes_.equal_range(dip);
			bool found = false;

			for (auto i = range.first; i != range.second; ++i)
			{
				if (i->second.second == dipole_id)
				{
					// there is already a matrix element for this dipole
					i->second.first.add(states.first);
					found = true;
				}
			}

			if (!found)
			{
				mes_.emplace(dip, std::make_pair(
					initial_state_set{states.first}, dipole_id));
			}

			// currently we only support one photon dipole per spectator
			if (photon_dipole_selected)
			{
				break;
			}
		}
		}
		}
		}
	}

	std::sort(dipoles_.begin(), dipoles_.end());
	dipoles_.shrink_to_fit();

	final_states_real_.shrink_to_fit();
	final_states_.shrink_to_fit();
	std::size_t const n = final_states_.size() + 2;
	ol_m2_.resize(n * (n - 1) / 2);
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

	auto const range = mes_.equal_range(dipole);
	auto const em = dipole.emitter();
	auto const sp = dipole.spectator();

	if (dipole.corr_type() == correction_type::qcd)
	{
		double as;
		ol.getparameter_double("alphas", &as);
		T const alphas = T(as);

		for (auto i = range.first; i != range.second; ++i)
		{
			auto const dipole_set = i->second.first;
			auto const id = i->second.second;
			auto const this_set = set.intersection(dipole_set);

			if (!this_set.empty())
			{
				double m2tree;
				double m2ew;
				ol.evaluate_cc(id, ol_phase_space_.data(), &m2tree,
					ol_m2_.data(), &m2ew);

				auto const k = std::min(em, sp);
				auto const l = std::max(em, sp);
				auto const index = k + l * (l - 1) / 2;

				T const result = alphas * T(ol_m2_.at(index));

				for (auto const state : this_set)
				{
					results.front().emplace_back(state, result);
				}
			}
		}
	}
	else if (dipole.corr_type() == correction_type::ew)
	{
		double a;
		ol.getparameter_double("alpha", &a);
		T const alpha = T(a);

		for (auto i = range.first; i != range.second; ++i)
		{
			auto const dipole_set = i->second.first;
			auto const id = i->second.second;
			auto const this_set = set.intersection(dipole_set);

			if (!this_set.empty())
			{
				double m2tree;
				ol.evaluate_tree(id, ol_phase_space_.data(), &m2tree);

				T const charge_em = charge_table_.at(id).at(em);
				T const charge_sp = charge_table_.at(id).at(sp);
				T const result = charge_em * charge_sp * alpha * T(m2tree);

				for (auto const state : this_set)
				{
					results.front().emplace_back(state, result);
				}
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
void ol_real_matrix_elements<T>::dipole_sc(
	hep::dipole const& dipole,
	std::vector<T> const& phase_space,
	std::array<T, 4> const& vector,
	hep::initial_state_set set,
	std::vector<hep::scales<T>> const& scales,
	std::vector<hep::initial_state_map<T>>& results_one,
	std::vector<hep::initial_state_map<T>>& results_two
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

	std::array<double, 4> double_vector = {
		static_cast <double> (vector.at(0)),
		static_cast <double> (vector.at(1)),
		static_cast <double> (vector.at(2)),
		static_cast <double> (vector.at(3))
	};

	auto const range = mes_.equal_range(dipole);
	auto const em = dipole.emitter();
	auto const sp = dipole.spectator();

	if (dipole.corr_type() == correction_type::qcd)
	{
		double as;
		ol.getparameter_double("alphas", &as);
		T const alphas = T(as);

		for (auto i = range.first; i != range.second; ++i)
		{
			auto const dipole_set = i->second.first;
			auto const id = i->second.second;
			auto const this_set = set.intersection(dipole_set);

			if (!this_set.empty())
			{
				double m2tree;
				double m2ew;
				ol.evaluate_cc(id, ol_phase_space_.data(), &m2tree,
					ol_m2_.data(), &m2ew);

				auto const k = std::min(em, sp);
				auto const l = std::max(em, sp);
				auto const index = k + l * (l - 1) / 2;

				T const result_one = alphas * T(ol_m2_.at(index));

				for (auto const state : this_set)
				{
					results_one.front().emplace_back(state, result_one);
				}

				ol.evaluate_sc(id, ol_phase_space_.data(), em + 1,
					double_vector.data(), ol_m2_.data());

				T const result_two = alphas * T(ol_m2_.at(sp));

				for (auto const state : this_set)
				{
					results_two.front().emplace_back(state, result_two);
				}
			}
		}
	}
	else if (dipole.corr_type() == correction_type::ew)
	{
		double a;
		ol.getparameter_double("alpha", &a);
		T const alpha = T(a);

		for (auto i = range.first; i != range.second; ++i)
		{
			auto const dipole_set = i->second.first;
			auto const id = i->second.second;
			auto const this_set = set.intersection(dipole_set);

			if (!this_set.empty())
			{
				T const charge_em = charge_table_.at(id).at(em);

				double m2tree;
				ol.evaluate_tree(id, ol_phase_space_.data(), &m2tree);

				T const result_one = -charge_em * charge_em * alpha * T(m2tree);

				for (auto const state : this_set)
				{
					results_one.front().emplace_back(state, result_one);
				}

				// OpenLoops returns the spin-correlator with an additonal sign
				ol.evaluate_sc(id, ol_phase_space_.data(), em + 1,
					double_vector.data(), ol_m2_.data());

				T const result_two = -charge_em * charge_em * alpha *
					T(-ol_m2_.at(sp));

				for (auto const state : this_set)
				{
					results_two.front().emplace_back(state, result_two);
				}
			}
		}
	}
	else
	{
		assert( false );
	}

	for (std::size_t i = 1; i != scales.size(); ++i)
	{
		results_one.at(i) = results_one.front();
		results_two.at(i) = results_two.front();
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

		// don't write zeros into `results`
		if (range.first == range.second)
		{
			continue;
		}

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
