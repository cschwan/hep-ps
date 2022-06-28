#include "hep/ps/correction_type.hpp"
#include "hep/ps/generate_dipole.hpp"
#include "hep/ps/me_type.hpp"
#include "hep/ps/ol_interface.hpp"
#include "hep/ps/ol_int_dipoles.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <map>
#include <stdexcept>
#include <utility>

namespace hep
{

template <typename T>
ol_int_dipoles<T>::ol_int_dipoles(
    std::vector<std::string> const& real_processes,
    std::vector<final_state> const& dipole_final_states,
    coupling_order order,
    dipole_veto const& veto,
    photon_dipole_selector const& selector
)
    : alphas_power_(order.alphas_power())
    , final_states_(dipole_final_states)
{
    auto& ol = ol_interface::instance();

    std::unordered_multimap<int_dipole, std::tuple<std::vector<int>, int, std::size_t, T>> mes;

    std::vector<int> dipole_ids;

    for (auto const& real_process : real_processes)
    {
        auto const& ids = ol_process_string_to_pdg_ids(real_process);

        // construct all possible EW and QCD dipoles
        for (auto const type : correction_type_list())
        {
        for (std::size_t i = 0; i != ids.size(); ++i)
        {
        for (std::size_t j = 2; j != ids.size(); ++j)
        {
            bool is_photon_dipole = false;
            bool photon_dipole_selected;

        for (std::size_t k = 0; k != ids.size(); ++k)
        {
            if ((i == j) || (i == k) || (j == k))
            {
                continue;
            }

            auto const type_i = pdg_id_to_particle_type(ids.at(i));
            auto const type_j = pdg_id_to_particle_type(ids.at(j));
            auto const type_k = pdg_id_to_particle_type(ids.at(k));
            auto const dip = dipole(i, j, k, type_i, type_j, type_k, type);

            auto splitting = generate_dipole(ids, dipole_ids, order, dip);

            if (dipole_ids.empty())
            {
                // there is no matrix element for this dipole
                continue;
            }

            if (veto(dipole_ids, final_states_, dip))
            {
                continue;
            }

            T const factor = T(final_state_symmetry_factor(dipole_ids)) /
                T(final_state_symmetry_factor(ids));

            photon_dipole_selected = false;

            if (type == correction_type::ew)
            {
                // if the unresolved particle is not a final-state photon, we only need one dipole
                if (ids.at(j) != pdg_id_of_photon())
                {
                    is_photon_dipole = true;

                    if (selector(dipole_ids, i, j, k))
                    {
                        photon_dipole_selected = true;
                    }
                    else
                    {
                        continue;
                    }

                    auto const photon_index = dip.emitter()
                        - ((dip.unresolved() < dip.emitter()) ? 1 : 0);

                    // OpenLoops quirk: a -> ff~ without alpha(0)/alpha by using `-22` for the
                    // photon
                    if (dipole_ids.at(photon_index) == pdg_id_of_photon())
                    {
                        dipole_ids.at(photon_index) = -pdg_id_of_photon();
                    }
                }
            }

            auto const ol_type = (type == correction_type::qcd) ?
                me_type::color_correlated : me_type::born;

            auto const process = pdg_ids_to_ol_process_string(dipole_ids);
            int const order_ew = (type == correction_type::qcd) ? order.alpha_power() :
                (order.alpha_power() - 1);
            int const order_qcd = (type == correction_type::qcd) ? (order.alphas_power() - 1) :
                order.alphas_power();
            int const dipole_id = ol.register_process(
                process.c_str(),
                ol_type,
                (ol_type == me_type::color_correlated) ? order.alphas_power() : order_qcd,
                (ol_type == me_type::color_correlated) ? order.alpha_power() : order_ew
            );

            std::size_t charge_table_index = -1;

            // add a charge table if neccessary
            if (type == correction_type::ew)
            {
                std::vector<T> charges;
                charges.reserve(dipole_ids.size());

                // we need the charges of the real matrix element
                for (int const id : dipole_ids)
                {
                    int const charge = pdg_id_to_charge_times_three(id);
                    charges.push_back(T(charge) / T(3.0));
                }

                // sign of the crossing
                charges.at(0) *= T(-1.0);
                charges.at(1) *= T(-1.0);

                // check if this charge table does already exist
                auto const result = std::find_if(charge_table_.begin(), charge_table_.end(),
                    [&](std::vector<T> const& table) {
                        return std::equal(table.begin(), table.end(), charges.begin());
                });

                if (result == charge_table_.end())
                {
                    charge_table_index = charge_table_.size();
                    charge_table_.emplace_back(std::move(charges));
                }
                else
                {
                    charge_table_index = std::distance(charge_table_.begin(), result);
                }
            }

            std::size_t const ins_i = (i < j) ? i : i - 1;
            std::size_t const ins_k = (k < j) ? k : k - 1;
            //auto const ins_type_i = pdg_id_to_particle_type(dipole_ids.at(ins_i));

            auto const term = int_dipole(ins_i, ins_k, splitting);

            // add a dipole if it doesn't exist yet
            if (std::find(dipoles_.begin(), dipoles_.end(), term) == dipoles_.end())
            {
                dipoles_.push_back(term);
            }

            mes.emplace(term, std::make_tuple(dipole_ids, dipole_id, charge_table_index, factor));

            if ((term.type() == insertion_term_type::initial_initial) ||
                (term.type() == insertion_term_type::initial_final))
            {
                auto const born = int_dipole{term.initial_particle(), splitting};

                if (std::find(dipoles_.begin(), dipoles_.end(), born) == dipoles_.end())
                {
                    dipoles_.push_back(born);
                }

                int const born_id = ol.register_process(process.c_str(), me_type::born,
                    order_qcd, order_ew);

                auto const range = mes.equal_range(born);
                auto const tuple = std::make_tuple(dipole_ids, born_id, charge_table_index, T(1.0));
                bool exists = false;

                for (auto i = range.first; i != range.second; ++i)
                {
                    if (std::get<1>(i->second) != std::get<1>(tuple) ||
                        (std::get<2>(i->second) != std::get<2>(tuple)))
                    {
                        continue;
                    }

                    if (!std::equal(std::get<0>(tuple).begin(), std::get<0>(tuple).end(),
                        std::get<0>(i->second).begin()))
                    {
                        continue;
                    }

                    // the term that we try to add already exists
                    exists = true;
                    break;
                }

                if (!exists)
                {
                    mes.emplace(born, tuple);
                }
            }

            // currently we only support one photon dipole per spectator
            if (photon_dipole_selected)
            {
                break;
            }
        }

            if (is_photon_dipole)
            {
                if (photon_dipole_selected)
                {
                    is_photon_dipole = false;
                }
                else
                {
                    throw std::runtime_error("dipole veto vetoes all photon dipoles");
                }
            }
        }
        }
        }
    }

    std::sort(dipoles_.begin(), dipoles_.end());
    dipoles_.shrink_to_fit();
    charge_table_.shrink_to_fit();
    mes_.reserve(mes.size());

    // we don't need the full particle info, only the initial state
    for (auto const& me : mes)
    {
        mes_.emplace(me.first, std::make_tuple(pdg_ids_to_states(std::get<0>(me.second)).first,
            std::get<1>(me.second), std::get<2>(me.second), std::get<3>(me.second)));
    }

    final_states_.shrink_to_fit();
    std::size_t const n = final_states_.size() + 2;
    ol_m2_.resize(n * (n - 1) / 2);
    ol_phase_space_.resize(5 * (n + 1));
}

template <typename T>
void ol_int_dipoles<T>::alphas(T alphas)
{
    auto& ol = ol_interface::instance();
    ol.setparameter_double("alphas", static_cast <double>(alphas));
}

template <typename T>
std::size_t ol_int_dipoles<T>::alphas_power() const
{
    return alphas_power_;
}

template <typename T>
void ol_int_dipoles<T>::correlated_me(
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
        ol_phase_space_.at(5 * i + 0) = static_cast <double> (phase_space.at(4 * i + 0));
        ol_phase_space_.at(5 * i + 1) = static_cast <double> (phase_space.at(4 * i + 1));
        ol_phase_space_.at(5 * i + 2) = static_cast <double> (phase_space.at(4 * i + 2));
        ol_phase_space_.at(5 * i + 3) = static_cast <double> (phase_space.at(4 * i + 3));
        ol_phase_space_.at(5 * i + 4) = 0.0;
    }

    std::size_t me = 0;

    for (auto const& dipole : dipoles_)
    {
        auto const range = mes_.equal_range(dipole);

        if (correction_type_of(dipole.vertex()) == correction_type::qcd)
        {
            double as;
            ol.getparameter_double("alphas", &as);
            T const alphas = T(as);

            for (auto i = range.first; i != range.second; ++i)
            {
                auto const state = std::get<0>(i->second);
                auto const id = std::get<1>(i->second);
                T const factor = std::get<3>(i->second);

                if (set.includes(state))
                {
                    if (id == -1)
                    {
                        results.at(me).emplace_back(state, T());
                        continue;
                    }

                    double m2tree;
                    double m2ew;
                    ol.evaluate_cc(id, ol_phase_space_.data(), &m2tree, ol_m2_, &m2ew);

                    if (dipole.type() == insertion_term_type::born)
                    {
                        T const casimir = casimir_operator<T>(state, dipole.initial_particle());
                        T const result = casimir * alphas * factor * T(m2tree);

                        results.at(me).emplace_back(state, result);
                    }
                    else
                    {
                        auto const em = dipole.emitter();
                        auto const sp = dipole.spectator();
                        auto const k = std::min(em, sp);
                        auto const l = std::max(em, sp);
                        auto const index = k + l * (l - 1) / 2;

                        T const result = alphas * factor * T(ol_m2_.at(index));

                        results.at(me).emplace_back(state, result);
                    }
                }
            }
        }
        else if (correction_type_of(dipole.vertex()) == correction_type::ew)
        {
            double a;
            ol.getparameter_double("alpha", &a);
            T const alpha = T(a);

            for (auto i = range.first; i != range.second; ++i)
            {
                auto const state = std::get<0>(i->second);
                auto const id = std::get<1>(i->second);
                auto const index = std::get<2>(i->second);
                T const factor = std::get<3>(i->second);

                if (set.includes(state))
                {
                    if (id == -1)
                    {
                        results.at(me).emplace_back(state, T());
                        continue;
                    }

                    double m2tree;
                    ol.evaluate_tree(id, ol_phase_space_.data(), &m2tree);

                    if (dipole.type() == insertion_term_type::born)
                    {
                        T charge;

                        if (dipole.vertex().unresolved() != pdg_id_of_photon())
                        {
                            charge = T(pdg_id_to_charge_times_three(dipole.vertex().unresolved()))
                                / T(3.0);
                        }
                        else
                        {
                            charge = charge_table_.at(index).at(dipole.initial_particle());
                        }

                        results.at(me).emplace_back(state, charge * charge * factor * T(alpha) * T(m2tree));
                    }
                    else
                    {
                        T charges;

                        if (dipole.vertex().unresolved() != pdg_id_of_photon())
                        {
                            charges = T(pdg_id_to_charge_times_three(dipole.vertex().unresolved()))
                                / T(3.0);
                            charges *= -charges;
                        }
                        else
                        {
                            auto const em = dipole.emitter();
                            auto const sp = dipole.spectator();
                            T const charge_em = charge_table_.at(index).at(em);
                            T const charge_sp = charge_table_.at(index).at(sp);

                            charges = charge_em * charge_sp;
                        }

                        T const result = charges * alpha * factor * T(m2tree);

                        results.at(me).emplace_back(state, result);
                    }
                }
            }
        }
        else
        {
            assert( false );
        }

        ++me;
    }
}

template <typename T>
std::vector<int_dipole> const& ol_int_dipoles<T>::dipoles() const
{
    return dipoles_;
}

template <typename T>
std::vector<final_state> const& ol_int_dipoles<T>::final_states() const
{
    return final_states_;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class ol_int_dipoles<double>;

}
