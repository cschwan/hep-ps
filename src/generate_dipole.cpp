#include "hep/ps/generate_dipole.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>
#include <cstdlib>

namespace hep
{

dipole_vertex generate_dipole(
    std::vector<int> const& process_pdg_ids,
    std::vector<int>& dipole_pdg_ids,
    coupling_order order,
    dipole const& dipole_info
) {
    auto const i = dipole_info.emitter();
    auto const j = dipole_info.unresolved();
    auto const k = dipole_info.spectator();
    auto const type = dipole_info.corr_type();

    int const id_i = process_pdg_ids.at(i);
    int const id_j = process_pdg_ids.at(j);
    int const id_k = process_pdg_ids.at(k);
    int const sign = ((i < 2) == (j < 2)) ? 1 : -1; // TODO: `j` should always be larger than 1

    dipole_vertex result;

    if (type == correction_type::ew)
    {
        bool const charged_i = pdg_id_has_charge(id_i);
        bool const charged_j = pdg_id_has_charge(id_j);
        bool const charged_k = pdg_id_has_charge(id_k);

        // fermion -> fermion + photon
        if (pdg_id_is_photon(id_j) && charged_i && charged_k)
        {
            dipole_pdg_ids = process_pdg_ids;
            dipole_pdg_ids.erase(dipole_pdg_ids.begin() + j);

            result = dipole_vertex(id_i, id_i, id_j);
        }
        // fermion -> photon + fermion
        else if (pdg_id_is_photon(id_i) && charged_j && (i < 2))
        {
            // TODO: NYI
            assert( false );
        }
        // photon -> fermion + antifermion
        else if (charged_j && ((id_i + sign * id_j) == 0) && ((i < 2) || (id_i > 0)))
        {
            dipole_pdg_ids = process_pdg_ids;
            dipole_pdg_ids.at(i) = pdg_id_of_photon();
            dipole_pdg_ids.erase(dipole_pdg_ids.begin() + j);

            result = dipole_vertex(pdg_id_of_photon(), id_i, pdg_id_anti(id_i));
        }
        else
        {
            dipole_pdg_ids.clear();
        }
    }
    else if (type == correction_type::qcd)
    {
        if (pdg_id_has_color(id_i) && pdg_id_has_color(id_j) && pdg_id_has_color(id_k))
        {
            // quark -> quark + gluon or gluon -> gluon + gluon
            if (pdg_id_is_gluon(id_j))
            {
                dipole_pdg_ids = process_pdg_ids;
                dipole_pdg_ids.erase(dipole_pdg_ids.begin() + j);

                result = dipole_vertex(id_i, id_i, id_j);
            }
            // quark -> gluon + quark
            else if (pdg_id_is_gluon(id_i) && (i < 2))
            {
                dipole_pdg_ids = process_pdg_ids;
                dipole_pdg_ids.at(i) = process_pdg_ids.at(j) * sign;
                dipole_pdg_ids.erase(dipole_pdg_ids.begin() + j);

                result = dipole_vertex(pdg_id_of_gluon(), id_i, pdg_id_anti(id_i));
            }
            // gluon -> quark + antiquark
            else if ((((id_i + sign * id_j) == 0) && ((i < 2) || (id_i > 0))))
            {
                dipole_pdg_ids = process_pdg_ids;
                dipole_pdg_ids.at(i) = pdg_id_of_gluon();
                dipole_pdg_ids.erase(dipole_pdg_ids.begin() + j);

                result = dipole_vertex(pdg_id_of_gluon(), id_i, pdg_id_anti(id_i));
            }
            else
            {
                dipole_pdg_ids.clear();
            }
        }
        else
        {
            dipole_pdg_ids.clear();
        }
    }
    else
    {
        assert( false );
    }

    std::size_t const gluons = std::count_if(dipole_pdg_ids.begin(), dipole_pdg_ids.end(),
        pdg_id_is_gluon);
    std::size_t const quarks = std::count_if(dipole_pdg_ids.begin(), dipole_pdg_ids.end(),
        pdg_id_is_quark);

    std::size_t n_qcd_min;
    std::size_t n_qcd_max;

    if (gluons >= (dipole_pdg_ids.size() - 2))
    {
        // pure gluon case and one quark line
        n_qcd_min = (gluons == dipole_pdg_ids.size()) ? (gluons - 2) : gluons;
        n_qcd_max = n_qcd_min;
    }
    else
    {
        // starting with two quark lines it is possible to have photons or
        // gluons in between the quark lines
        n_qcd_min = gluons;
        n_qcd_max = n_qcd_min + quarks / 2;
    }

    std::size_t const order_qcd = (type == correction_type::qcd) ?
        (order.alphas_power() - 1) : order.alphas_power();

    // check if the matrix element with given coupling order exists
    if ((order_qcd < n_qcd_min) || (order_qcd > n_qcd_max))
    {
        dipole_pdg_ids.clear();
    }

    return result;
}

}
