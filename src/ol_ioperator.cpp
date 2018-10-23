#include "hep/ps/ol_ioperator.hpp"
#include "hep/ps/ol_interface.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <stdexcept>

namespace hep
{

template <typename T>
ol_ioperator<T>::ol_ioperator(
    std::string const& process,
    coupling_order const& order,
    regularization_scheme scheme
)
    : order_{order}
{
    auto& ol = ol_interface::instance();

    // TODO: fix the determination of the coupling order
    ol.setparameter_int("order_ew", order.alpha_power());
//    ol.setparameter_int("order_qcd", order.alphas_power());
//    ol.setparameter_int("verbose", 3);

    // evaluates the I-operator
    ol.setparameter_int("IR_on", 2);
    // we don't need the stability system for the I-operator
    ol.setparameter_int("stability_mode", 11);

    switch (scheme)
    {
    case regularization_scheme::dim_reg_blha:
        ol.setparameter_int("polenorm", 0);
        break;

    case regularization_scheme::dim_reg_coli:
        ol.setparameter_int("polenorm", 1);
        break;

    default:
        assert( false );
    }

    id_ = ol.register_process(process.c_str(), 11);

    if (id_ == -1)
    {
        throw std::runtime_error("OpenLoops didn't find the given process");
    }

    auto const& pdg_ids = ol_process_string_to_pdg_ids(process);
    state_ = partons_to_initial_state(pdg_id_to_parton(pdg_ids.at(0)),
        pdg_id_to_parton(pdg_ids.at(1)));
    final_states_.reserve(pdg_ids.size() - 2);

    for (std::size_t i = 2; i != pdg_ids.size(); ++i)
    {
        final_states_.push_back(pdg_id_to_final_state(pdg_ids.at(i)));
    }
}

template <typename T>
void ol_ioperator<T>::alphas(T alphas)
{
    auto& ol = ol_interface::instance();

    ol.setparameter_double("alphas", static_cast <double> (alphas));
}

template <typename T>
std::size_t ol_ioperator<T>::alphas_power() const
{
    return order_.alphas_power();
}

template <typename T>
void ol_ioperator<T>::borns(
    std::vector<T> const& phase_space,
    initial_state_set /*set*/,
    nonstd::span<scales<T> const> scales,
    std::vector<initial_state_map<T>>& results
) {
    auto& ol = ol_interface::instance();

    // TODO: set also renormalization scale and loop over it
    ol.setparameter_double("mu", static_cast <double> (scales[0].regularization()));

    std::size_t n = ol.n_external(id_);
    std::vector<double> pp(5 * n);

    for (std::size_t i = 0; i != n; ++i)
    {
        pp.at(5 * i + 0) = static_cast <double> (phase_space.at(4 * i + 0));
        pp.at(5 * i + 1) = static_cast <double> (phase_space.at(4 * i + 1));
        pp.at(5 * i + 2) = static_cast <double> (phase_space.at(4 * i + 2));
        pp.at(5 * i + 3) = static_cast <double> (phase_space.at(4 * i + 3));

        // supposed to be the mass, is currently ignored by OpenLoops
        pp.at(5 * i + 4) = 0.0;
    }

    double m2tree;
    double m2loop[3];
    double m2ir1[3];
    double m2loop2[5];
    double m2ir2[5];
    double acc;

    ol.evaluate_full(id_, pp.data(), &m2tree, m2loop, m2ir1, m2loop2, m2ir2, &acc);

    for (auto& result : results)
    {
        result.emplace_back(state_, T(m2ir1[0]));
    }
}

template <typename T>
std::vector<hep::final_state> const& ol_ioperator<T>::final_states() const
{
    return final_states_;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class ol_ioperator<double>;

}
