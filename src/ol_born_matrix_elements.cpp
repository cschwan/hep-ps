#include "hep/ps/ol_born_matrix_elements.hpp"
#include "hep/ps/ol_interface.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>
#include <stdexcept>

namespace hep
{

template <typename T>
ol_born_matrix_elements<T>::ol_born_matrix_elements(
	std::vector<std::string> const& processes,
	std::size_t alphas_power
)
	: alphas_power_(alphas_power)
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

		ids_.emplace(state, ol.register_process(process.c_str(), 1));
	}

	final_states_.shrink_to_fit();
	ol_phase_space_.resize(5 * (final_states_.size() + 2));
}

template <typename T>
void ol_born_matrix_elements<T>::alphas(T alphas)
{
	auto& ol = ol_interface::instance();
	ol.setparameter_double("alphas", static_cast <double>(alphas));
}

template <typename T>
std::size_t ol_born_matrix_elements<T>::alphas_power() const
{
	return alphas_power_;
}

template <typename T>
void ol_born_matrix_elements<T>::borns(
	std::vector<T> const& phase_space,
	hep::initial_state_set set,
	std::vector<scales<T>> const& scales,
	std::vector<initial_state_map<T>>& results
) {
	auto& ol = hep::ol_interface::instance();
	auto& result = results.front();

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
		auto const range = ids_.equal_range(state);
		T value = T();

		for (auto i = range.first; i != range.second; ++i)
		{
			ol.evaluate_tree(i->second, ol_phase_space_.data(), &m2tree);
			value += T(m2tree);
		}

		result.emplace_back(state, value);
	}

	for (std::size_t i = 1; i != scales.size(); ++i)
	{
		results.at(i) = result;
	}
}

template <typename T>
std::vector<final_state> const& ol_born_matrix_elements<T>::final_states() const
{
	return final_states_;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class ol_born_matrix_elements<double>;

}
