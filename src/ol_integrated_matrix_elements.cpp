#include "hep/ps/ol_integrated_matrix_elements.hpp"
#include "hep/ps/ol_interface.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace hep
{

template <typename T>
ol_integrated_matrix_elements<T>::ol_integrated_matrix_elements(
	std::vector<std::string> const& processes,
	std::size_t alphas_power,
	correction_type type
)
	: alphas_power_(alphas_power)
	, type_(type)
{
	auto& ol = ol_interface::instance();
	ol.setparameter_int("order_qcd", alphas_power);

	for (auto const& process : processes)
	{
		auto const& pdg_ids = ol_process_string_to_pdg_ids(process);

		std::vector<final_state> final_states(pdg_ids.size() - 2);
		std::transform(pdg_ids.begin() + 2, pdg_ids.end(), final_states.begin(),
			pdg_id_to_final_state);

		auto const state = partons_to_initial_state(
			pdg_id_to_parton(pdg_ids.at(0)), pdg_id_to_parton(pdg_ids.at(1)));

		if (final_states_.empty())
		{
			final_states_ = final_states;
		}
		else
		{
			if (!std::equal(final_states_.begin(), final_states_.end(),
				final_states.begin()))
			{
				throw std::invalid_argument("processes are not compatible");
			}
		}

		int const process_id = ol.register_process(process.c_str(), 1);
		ids_.emplace(state, process_id);

		if (type == correction_type::qcd)
		{
			// TODO: NYI
			assert( false );
		}
		else if (type == correction_type::ew)
		{
			std::vector<T> charges(pdg_ids.size());
			std::transform(pdg_ids.begin(), pdg_ids.end(), charges.begin(),
				pdg_id_to_charge<T>);

			charges.at(0) *= T(-1.0);
			charges.at(1) *= T(-1.0);

			std::vector<std::size_t> indices;
			indices.reserve(charges.size());

			for (std::size_t i = 0; i != charges.size(); ++i)
			{
				if (charges.at(i) != T())
				{
					indices.push_back(i);
				}
			}

			charge_table_.emplace(process_id, std::move(charges));

			terms_.emplace_back(0);
			terms_.emplace_back(1);

			for (std::size_t i = 0; i < indices.size() - 1; ++i)
			{
				for (std::size_t j = i + 1; j != indices.size(); ++j)
				{
					auto const k = indices.at(i);
					auto const l = indices.at(j);
					auto const type_k = pdg_id_to_particle_type(pdg_ids.at(k));
					auto const type_l = pdg_id_to_particle_type(pdg_ids.at(l));

					terms_.emplace_back(k, type_k, l);
					terms_.emplace_back(l, type_l, k);
				}
			}
		}
	}

	std::sort(terms_.begin(), terms_.end());
	terms_.erase(std::unique(terms_.begin(), terms_.end()), terms_.end());
	terms_.shrink_to_fit();

	final_states_.shrink_to_fit();
	ol_phase_space_.resize(5 * (final_states_.size() + 2));
}

template <typename T>
void ol_integrated_matrix_elements<T>::alphas(T alphas)
{
	auto& ol = ol_interface::instance();
	ol.setparameter_double("alphas", static_cast <double>(alphas));
}

template <typename T>
std::size_t ol_integrated_matrix_elements<T>::alphas_power() const
{
	return alphas_power_;
}

template <typename T>
void ol_integrated_matrix_elements<T>::correlated_me(
	std::vector<T> const& phase_space,
	initial_state_set set,
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
		// TODO: NYI
		assert( false );
	}
	else if (type_ == correction_type::ew)
	{
		double m2tree;
		double alpha;
		ol.getparameter_double("alpha", &alpha);

		for (auto const state : set)
		{
			auto const range = ids_.equal_range(state);

			for (auto i = range.first; i != range.second; ++i)
			{
				ol.evaluate_tree(i->second, ol_phase_space_.data(), &m2tree);

				for (std::size_t j = 0; j != terms_.size(); ++j)
				{
					auto const& term = terms_.at(j);

					if (term.type() == insertion_term_type::born)
					{
						T const charge = charge_table_.at(i->second).at(
							term.initial_particle());

						results.at(j).emplace_back(state,
							charge * charge * T(alpha) * T(m2tree));
					}
					else
					{
						T const charge_em = charge_table_.at(i->second).at(
							term.emitter());
						T const charge_sp = charge_table_.at(i->second).at(
							term.spectator());

						results.at(j).emplace_back(state,
							charge_em * charge_sp * T(alpha) * T(m2tree));
					}
				}
			}
		}
	}
	else
	{
		assert( false );
	}
}

template <typename T>
std::vector<final_state> const&
ol_integrated_matrix_elements<T>::final_states() const
{
	return final_states_;
}

template <typename T>
std::vector<insertion_term> const&
ol_integrated_matrix_elements<T>::insertion_terms() const
{
	return terms_;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class ol_integrated_matrix_elements<double>;

}
