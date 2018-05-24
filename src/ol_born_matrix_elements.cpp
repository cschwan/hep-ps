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
	std::size_t alphas_power,
	bool loop_mes,
	regularization_scheme scheme,
	bool set_ol_coupling_order
)
	: alphas_power_(alphas_power)
	, loop_mes_(loop_mes)
{
	auto& ol = ol_interface::instance();

	if (set_ol_coupling_order)
	{
		ol.setparameter_int("order_qcd", alphas_power);
	}

	if (loop_mes)
	{
		// TODO: if we calculate EW corrections, set `ew_renorm` to `1`

		// we will calculate the counterterms separately
		ol.setparameter_int("CT_on", 0);

		if (scheme == regularization_scheme::dim_reg_blha)
		{
			ol.setparameter_int("polenorm", 0);
		}
		else if (scheme == regularization_scheme::dim_reg_coli)
		{
			ol.setparameter_int("polenorm", 1);
		}
		else
		{
			assert( false );
		}
	}

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

		int const amptype = loop_mes ? 11 : 1;
		ids_.emplace(state, ol.register_process(process.c_str(), amptype));
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

	if (loop_mes_)
	{
		double m2tree;
		double m2loop[3];
		double acc;

		double const mureg =
			static_cast <double> (scales.front().regularization());

		ol.setparameter_double("mureg", mureg);
		ol.setparameter_double("muren",
			static_cast <double> (scales.front().renormalization()));

		// TODO: `scales` usually has more than one entry with the same
		// renormalization scale - cache it?

		for (auto const state : set)
		{
			auto const range = ids_.equal_range(state);
			T value = T();

			if (range.first == range.second)
			{
				continue;
			}

			for (auto i = range.first; i != range.second; ++i)
			{
				// TODO: this assumes that the matrix element itself does
				// not depend directly on the renormalization scale, for
				// example through MSBar masses

				ol.evaluate_loop(i->second, ol_phase_space_.data(), &m2tree,
					m2loop, &acc);
				value += T(m2loop[0]);
			}

			for (std::size_t i = 0; i != scales.size(); ++i)
			{
				results.at(i).emplace_back(state, value);
			}
		}

		for (std::size_t j = 0; j != scales.size(); ++j)
		{
			if (scales.at(j).regularization() !=
				scales.front().regularization())
			{
				throw std::runtime_error(
					"regularization scales must be the same");
			}

			double const muren =
				static_cast <double> (scales.at(j).renormalization());

			ol.setparameter_double("muren", muren);

			std::size_t k = 0;

			for (auto const state : set)
			{
				auto const range = ids_.equal_range(state);
				T value = T();

				for (auto i = range.first; i != range.second; ++i)
				{
					ol.evaluate_ct(i->second, ol_phase_space_.data(), &m2tree,
						m2loop);
					value += T(m2loop[0]);
				}

				results.at(j).at(k++).second += value;
			}
		}
	}
	else
	{
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
			results.at(i) = results.front();
		}
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
