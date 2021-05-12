#include "hep/ps/me_type.hpp"
#include "hep/ps/ol_born_matrix_elements.hpp"
#include "hep/ps/ol_interface.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace hep
{

template <typename T>
ol_born_matrix_elements<T>::ol_born_matrix_elements(
    std::vector<std::string> const& processes,
    coupling_order const& order,
    bool loop_mes,
    regularization_scheme scheme
)
    : alphas_power_{order.alphas_power()}
    , loop_mes_(loop_mes)
{
    auto& ol = ol_interface::instance();

    if (loop_mes)
    {
        // TODO: if we calculate EW corrections, set `ew_renorm` to `1`

        // we will calculate the counterterms separately
        //ol.setparameter_int("CT_on", 0);

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
        auto const& states = pdg_ids_to_states(pdg_ids);

        if (final_states_.empty())
        {
            final_states_ = states.second;
        }
        else
        {
            if (!std::equal(states.second.begin(), states.second.end(),
                final_states_.begin()))
            {
                throw std::invalid_argument("processes are not compatible");
            }
        }

        auto amptype = loop_mes ? me_type::loop : me_type::born;
        ids_.emplace(states.first, ol.register_process(process.c_str(), amptype,
            order.alphas_power(), order.alpha_power()));
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
    nonstd::span<hep::scales<T> const> scales,
    std::vector<initial_state_map<T>>& results
) {
    auto& ol = hep::ol_interface::instance();
    auto& result = results.front();

    // TODO: for the time being we assume that the matrix elements do only
    // indirectly depend on the renormalization scales (through alphas)

    std::size_t const n = phase_space.size() / 4;

    for (std::size_t i = 0; i != n; ++i)
    {
        ol_phase_space_.at(5 * i + 0) = static_cast <double> (phase_space.at(4 * i + 0));
        ol_phase_space_.at(5 * i + 1) = static_cast <double> (phase_space.at(4 * i + 1));
        ol_phase_space_.at(5 * i + 2) = static_cast <double> (phase_space.at(4 * i + 2));
        ol_phase_space_.at(5 * i + 3) = static_cast <double> (phase_space.at(4 * i + 3));
        ol_phase_space_.at(5 * i + 4) = 0.0;
    }

    if (loop_mes_)
    {
        double m2tree;
        double m2loop[3];
        double acc;

        double const mureg = static_cast <double> (scales[0].regularization());

        ol.setparameter_double("mureg", mureg);
        ol.setparameter_double("muren", static_cast <double> (scales[0].renormalization()));

        // TODO: `scales` usually has more than one entry with the same renormalization scale -
        // cache it?

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
                // TODO: this assumes that the matrix element itself does not depend directly on the
                // renormalization scale, for example through MSBar masses

                ol.evaluate_loop(i->second, ol_phase_space_.data(), &m2tree, m2loop, &acc);

                double const k = std::fabs(m2loop[0] / m2tree);

                if (k > 1000.0)
                {
                    auto const flags = std::cout.flags();
                    auto const precision = std::cout.precision();

                    std::cout << "virtual PS cut (K=" << k << "): ";
                    std::cout << std::scientific << std::setprecision(17)
                        << ol_phase_space_[ 0] << ", " << ol_phase_space_[ 1] << ", " << ol_phase_space_[ 2] << ", " << ol_phase_space_[ 3] << ", "
                        << ol_phase_space_[ 4] << ", " << ol_phase_space_[ 5] << ", " << ol_phase_space_[ 6] << ", " << ol_phase_space_[ 7] << ", "
                        << ol_phase_space_[ 8] << ", " << ol_phase_space_[ 9] << ", " << ol_phase_space_[10] << ", " << ol_phase_space_[11] << ", "
                        << ol_phase_space_[12] << ", " << ol_phase_space_[13] << ", " << ol_phase_space_[14] << ", " << ol_phase_space_[15] << ", "
                        << ol_phase_space_[16] << ", " << ol_phase_space_[17] << ", " << ol_phase_space_[18] << ", " << ol_phase_space_[19] << ", "
                        << ol_phase_space_[20] << ", " << ol_phase_space_[21] << ", " << ol_phase_space_[22] << ", " << ol_phase_space_[23] << ", "
                        << ol_phase_space_[24] << ", " << ol_phase_space_[25] << ", " << ol_phase_space_[26] << ", " << ol_phase_space_[27] << ", "
                        << ol_phase_space_[28] << ", " << ol_phase_space_[29] << ", " << ol_phase_space_[30] << ", " << ol_phase_space_[31] << std::endl;

                    std::cout.flags(flags);
                    std::cout.precision(precision);

                    // cut the entire PS point
                    for (std::size_t i = 0; i != size(scales); ++i)
                    {
                        results.at(i).clear();
                    }

                    return;
                }

                value += T(m2loop[0]);
            }

            for (std::size_t i = 0; i != size(scales); ++i)
            {
                results.at(i).emplace_back(state, value);
            }
        }

        for (std::size_t j = 0; j != size(scales); ++j)
        {
            if (scales[j].regularization() != scales[0].regularization())
            {
                throw std::runtime_error("regularization scales must be the same");
            }

            double muren = 0.0;

            std::size_t k = 0;

            for (auto const state : set)
            {
                auto const range = ids_.equal_range(state);

                if (range.first == range.second)
                {
                    continue;
                }

                T value = T();

                for (auto i = range.first; i != range.second; ++i)
                {
                    muren = static_cast <double> (scales[j].renormalization());
                    ol.setparameter_double("muren", muren);
                    ol.evaluate_ct(i->second, ol_phase_space_.data(), &m2tree, m2loop);
                    value += T(m2loop[0]);

                    double m2loop0[3];

                    muren = static_cast <double> (scales[0].renormalization());
                    ol.setparameter_double("muren", muren);
                    ol.evaluate_ct(i->second, ol_phase_space_.data(), &m2tree, m2loop0);
                    value -= T(m2loop0[0]);
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

            if (range.first == range.second)
            {
                continue;
            }

            T value = T();

            for (auto i = range.first; i != range.second; ++i)
            {
                ol.evaluate_tree(i->second, ol_phase_space_.data(), &m2tree);
                value += T(m2tree);
            }

            result.emplace_back(state, value);
        }

        for (std::size_t i = 1; i != size(scales); ++i)
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
