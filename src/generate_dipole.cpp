#include "hep/ps/generate_dipole.hpp"
#include "hep/ps/pdg_functions.hpp"

#include <algorithm>

namespace hep
{

std::vector<int> generate_dipole(
	std::vector<int> const& process_pdg_ids,
	hep::coupling_order order,
	dipole const& dipole_info
) {
	auto const i = dipole_info.emitter();
	auto const j = dipole_info.unresolved();
	auto const k = dipole_info.spectator();
	auto const type = dipole_info.corr_type();

	std::vector<int> result;

	int const id_i = process_pdg_ids.at(i);
	int const id_j = process_pdg_ids.at(j);
	int const id_k = process_pdg_ids.at(k);
	int const sign = ((i < 2) == (j < 2)) ? 1 : -1;

	if (type == hep::correction_type::ew)
	{
		bool const charged_i = hep::pdg_id_has_charge(id_i);
		bool const charged_j = hep::pdg_id_has_charge(id_j);
		bool const charged_k = hep::pdg_id_has_charge(id_k);

		// fermion -> fermion + photon
		if (hep::pdg_id_is_photon(id_j) && charged_i && charged_k)
		{
			result = process_pdg_ids;
			result.erase(result.begin() + j);
		}
		// fermion -> photon + fermion
		else if (hep::pdg_id_is_photon(id_i) && charged_j && (i < 2))
		{
			// TODO: NYI
			assert( false );
		}
		// photon -> fermion + antifermion
		else if (charged_j && ((id_i + sign * id_j) == 0))
		{
			result = process_pdg_ids;
			result.at(i) = hep::pdg_id_of_photon();
			result.erase(result.begin() + j);
		}
	}
	else if (type == hep::correction_type::qcd)
	{
		if (hep::pdg_id_has_color(id_i) && hep::pdg_id_has_color(id_j) &&
			hep::pdg_id_has_color(id_k))
		{
			// quark -> quark + gluon or gluon -> gluon + gluon
			if (hep::pdg_id_is_gluon(id_j))
			{
				result = process_pdg_ids;
				result.erase(result.begin() + j);
			}
			// fermion -> gluon + fermion
			else if (hep::pdg_id_is_gluon(id_i) && (i < 2))
			{
				result = process_pdg_ids;
				result.at(i) = process_pdg_ids.at(j) * sign;
				result.erase(result.begin() + j);
			}
			// gluon -> fermion + antifermion
			else if ((id_i + sign * id_j) == 0)
			{
				result = process_pdg_ids;
				result.at(i) = hep::pdg_id_of_gluon();
				result.erase(result.begin() + j);
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
		(order.alphas_power() - 1) : order.alphas_power();

	// check if the matrix element with given coupling order exists
	if ((order_qcd < n_qcd_min) || (order_qcd > n_qcd_max))
	{
		result.clear();
	}

	return result;
}

}
