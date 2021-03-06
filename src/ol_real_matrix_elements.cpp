#include "hep/ps/correction_type.hpp"
#include "hep/ps/generate_dipole.hpp"
#include "hep/ps/me_type.hpp"
#include "hep/ps/ol_interface.hpp"
#include "hep/ps/ol_real_matrix_elements.hpp"
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
ol_real_matrix_elements<T>::ol_real_matrix_elements(
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

    std::vector<int> dipole_ids;

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

        int const real_id = ol.register_process(process.c_str(), me_type::born,
            order.alphas_power(), order.alpha_power());
        ids_reals_.emplace(states.first, real_id);

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

            auto const splitting = generate_dipole(ids, dipole_ids, order, dip);

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

            me_type ol_type = me_type::born;

            if (pdg_id_to_particle_type(splitting.internal()) == particle_type::boson)
            {
                // spin-correlated ME
                ol_type = me_type::spin_correlated;
            }
            else if (type == correction_type::qcd)
            {
                // color-correlated ME
                ol_type = me_type::color_correlated;
            }

            // if the dipole is a charge-correlator, we construct it from a Born matrix element; in
            // that case we have to modify the coupling orders; if it's a colour correlator instead
            // we get it from a loop matrix element for which we use the coupling order of the real
            // matrix element

            auto const process = pdg_ids_to_ol_process_string(dipole_ids);
            int const order_ew = (type == correction_type::qcd) ? order.alpha_power() :
                (order.alpha_power() - 1);
            int const order_qcd = (type == correction_type::qcd) ? order.alphas_power() :
                order.alphas_power();
            int const dipole_id = ol.register_process(process.c_str(), ol_type, order_qcd,
                order_ew);

            int charge_table_index = -1;

            // add a charge table if neccessary
            if (type == correction_type::ew)
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

            // add a dipole if it doesn't exist yet
            if (std::find(dipoles_.begin(), dipoles_.end(), dip) == dipoles_.end())
            {
                dipoles_.push_back(dip);
            }

            auto range = mes_.equal_range(dip);
            bool found = false;

            for (auto i = range.first; i != range.second; ++i)
            {
                // check if there is already a dipole that has the same ME
                if ((std::get<0>(i->second) == states.first) &&
                    (std::get<1>(i->second) == dipole_id) &&
                    (std::get<3>(i->second) == charge_table_index))
                {
                    std::get<2>(i->second) += factor;

                    found = true;

                    break;
                }
            }

            if (!found)
            {
                mes_.emplace(dip, std::make_tuple(states.first, dipole_id, factor,
                    charge_table_index));
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
    nonstd::span<scales<T> const> scales,
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

    auto const range = mes_.equal_range(dipole);
    auto const un = dipole.unresolved();
    auto const em_real = dipole.emitter();
    auto const sp_real = dipole.spectator();

    // when we remove the particle with index `un`, the emitter and spectator
    // indices for the born change
    auto const em_born = em_real - (em_real > un ? 1 : 0);
    auto const sp_born = sp_real - (sp_real > un ? 1 : 0);

    if (dipole.corr_type() == correction_type::qcd)
    {
        double as;
        ol.getparameter_double("alphas", &as);
        T const alphas = T(as);

        for (auto i = range.first; i != range.second; ++i)
        {
            auto const state = std::get<0>(i->second);
            auto const id = std::get<1>(i->second);
            auto const final_state_factor = std::get<2>(i->second);

            if (set.includes(state))
            {
                double m2tree;
                double m2ew;
                ol.evaluate_cc(id, ol_phase_space_.data(), &m2tree, ol_m2_.data(), &m2ew);

                auto const k = std::min(em_born, sp_born);
                auto const l = std::max(em_born, sp_born);
                auto const index = k + l * (l - 1) / 2;

                T const result = final_state_factor * alphas * T(ol_m2_.at(index));

                results.front().emplace_back(state, result);
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
            auto const state = std::get<0>(i->second);
            auto const id = std::get<1>(i->second);
            auto const final_state_factor = std::get<2>(i->second);
            auto const index = std::get<3>(i->second);

            if (set.includes(state))
            {
                double m2tree;
                ol.evaluate_tree(id, ol_phase_space_.data(), &m2tree);

                T const charge_un = charge_table_.at(index).at(dipole.unresolved());
                T const charge_em = charge_table_.at(index).at(em_real);

                T charges;

                if (charge_un != T())
                {
                    charges = -charge_un * charge_un;

                    if (charge_em == T())
                    {
                        // multiply with the number of colors
                        charges *= T(3.0);
                    }
                }
                else
                {
                    T const charge_sp = charge_table_.at(index).at(sp_real);

                    charges = charge_em * charge_sp;
                }

                T const result = final_state_factor * charges * alpha * T(m2tree);

                results.front().emplace_back(state, result);
            }
        }
    }
    else
    {
        assert( false );
    }

    for (std::size_t i = 1; i != size(scales); ++i)
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
    nonstd::span<hep::scales<T> const> scales,
    std::vector<hep::initial_state_map<T>>& results_one,
    std::vector<hep::initial_state_map<T>>& results_two
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

    std::array<double, 4> double_vector = {
        static_cast <double> (vector.at(0)),
        static_cast <double> (vector.at(1)),
        static_cast <double> (vector.at(2)),
        static_cast <double> (vector.at(3))
    };

    auto const range = mes_.equal_range(dipole);
    auto const un = dipole.unresolved();
    auto const em_real = dipole.emitter();
    auto const sp_real = dipole.spectator();

    // when we remove the particle with index `un`, the emitter and spectator
    // indices for the born change
    auto const em_born = em_real - (em_real > un ? 1 : 0);
    auto const sp_born = sp_real - (sp_real > un ? 1 : 0);

    if (dipole.corr_type() == correction_type::qcd)
    {
        double as;
        ol.getparameter_double("alphas", &as);
        T const alphas = T(as);

        for (auto i = range.first; i != range.second; ++i)
        {
            auto const state = std::get<0>(i->second);
            auto const id = std::get<1>(i->second);
            auto const final_state_factor = std::get<2>(i->second);

            if (set.includes(state))
            {
                double m2tree;
                double m2ew;
                ol.evaluate_cc(id, ol_phase_space_.data(), &m2tree, ol_m2_.data(), &m2ew);

                auto const k = std::min(em_born, sp_born);
                auto const l = std::max(em_born, sp_born);
                auto const index = k + l * (l - 1) / 2;

                T const result_one = final_state_factor * alphas * T(ol_m2_.at(index));

                results_one.front().emplace_back(state, result_one);

                ol.evaluate_sc(id, ol_phase_space_.data(), em_born + 1, double_vector.data(),
                    ol_m2_.data());

                T const result_two = final_state_factor * alphas * T(ol_m2_.at(sp_born));

                results_two.front().emplace_back(state, result_two);
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
            auto const state = std::get<0>(i->second);
            auto const id = std::get<1>(i->second);
            auto const final_state_factor = std::get<2>(i->second);
            auto const index = std::get<3>(i->second);

            if (set.includes(state))
            {
                T const charge_em = charge_table_.at(index).at(em_real);
                T const nc = T(3.0);

                // `em_real >= 2` and `em_born >= 2` are the same statements,
                // because `un` is always larger or equal to two
                T const initial_state_factor = ((em_real >= 2) ? nc : T(1.0)) *
                    charge_em * charge_em;

                double m2tree;
                ol.evaluate_tree(id, ol_phase_space_.data(), &m2tree);

                T const result_one = -final_state_factor * initial_state_factor * alpha * T(m2tree);

                results_one.front().emplace_back(state, result_one);

                // OpenLoops returns the spin-correlator with an additonal sign
                ol.evaluate_sc(id, ol_phase_space_.data(), em_born + 1, double_vector.data(),
                    ol_m2_.data());

                T const result_two = -final_state_factor * initial_state_factor * alpha *
                    T(-ol_m2_.at(sp_born));

                results_two.front().emplace_back(state, result_two);
            }
        }
    }
    else
    {
        assert( false );
    }

    for (std::size_t i = 1; i != size(scales); ++i)
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
    nonstd::span<scales<T> const> scales,
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

    for (std::size_t i = 1; i != size(scales); ++i)
    {
        results.at(i) = results.front();
    }
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class ol_real_matrix_elements<double>;

}
