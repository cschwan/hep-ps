#include "hep/ps/inv_idx.hpp"
#include "hep/ps/lusifer_diagram_generator.hpp"

#include "lusifer_interfaces.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>

namespace hep
{

template <typename T>
std::vector<diagram> lusifer_diagram_generator(
	std::vector<std::string> const& processes,
	lusifer_constants<T> const& constants
) {
	int maxex = -1;
	int maxgen = -1;
	lusifer_extra_max(&maxex, &maxgen);

	int nex = 0;

	// FORTRAN counting: use the first generator
	int g = 1;

	for (auto const& process : processes)
	{
		// each particle must be specified with three characters
		assert( process.size() % 3 == 0 );

		if (nex == 0)
		{
			// there must be at least four particles
			assert( process.size() >= 4 * 3 );

			nex = process.size() / 3;

			double mw = constants.mass_w;
			double gw = constants.width_w;
			double mz = constants.mass_z;
			double gz = constants.width_z;
			double mh = constants.mass_h;
			double gh = constants.width_h;
			double mt = constants.mass_t;
			double gt = constants.width_t;

			// set constants and the number of particles
			lusifer_extra_set(g, nex, mw, gw, mz, gz, mh, gh, mt, gt);

			// the number of particles must be supported by the generator
			assert( nex <= maxex );
		}

		// all processes must have the same number of particles
		assert( nex == static_cast <int> (process.size() / 3) );

		// fill up the string with three spaces for particle not used
		std::string process0 = process;
		process0.append(3 * (maxex - nex), ' ');

		// do not couple the Higgs to light fermions
		int lightfermions = 0;
		// do not include cuts in the phase space generation
		int includecuts = 0;
		// do not print channel information
		int sout = 0;

		lusifer_initphasespace(
			process0.c_str(),
			g,
			lightfermions,
			includecuts,
			sout
		);
	}

	int channels;
	lusifer_extra_data(g, &channels);

	// there must be at least one channel, otherwise something went wrong
	assert( channels > 0 );

	std::vector<diagram> diagrams;
	diagrams.reserve(channels);

	for (int i = 0; i != channels; ++i)
	{
		std::vector<vertex> vertices;
		std::vector<propagator> propagators;

		// the number of decays is always smaller or equal than the number of
		// invariants because the PS generator must generate invariants that do
		// not belong to a decay in the diagram
		int const decays = lusifer_cdecay.ndecay[0][i];

		for (int j = 0; j != decays; ++j)
		{
			vertices.emplace_back(
				inv_idx(lusifer_cdecay.indecay[0][i][j], nex),
				inv_idx(lusifer_cdecay.out1decay[0][i][j], nex),
				inv_idx(lusifer_cdecay.out2decay[0][i][j], nex)
			);

			auto const begin = std::begin(lusifer_cinv.ininv[0][i]);
			auto const end = begin + decays;

			// find the invariant corresponding to the decay
			auto const result = std::find(begin, end,
				lusifer_cdecay.indecay[0][i][j]);

			// there must be a valid search result, otherwise something went
			// wrong in the initialization of the phase space generator
			assert( result != end );

			propagators.emplace_back(
				*result,
				inv_idx(lusifer_cdecay.indecay[0][i][j], nex)
			);
		}

		for (int j = 0; j != lusifer_cprocess.nprocess[0][i]; ++j)
		{
			vertices.emplace_back(
				inv_idx(lusifer_cprocess.in1process[0][i][j], nex),
				inv_idx(lusifer_cprocess.out1process[0][i][j], nex),
				inv_idx(lusifer_cprocess.virtprocess[0][i][j], nex)
			);

			propagators.emplace_back(
				lusifer_cprocess.idhepprocess[0][i][j],
				inv_idx(lusifer_cprocess.virtprocess[0][i][j], nex)
			);
		}

		int last = lusifer_cprocess.nprocess[0][i] - 1;

		vertices.emplace_back(
			inv_idx(lusifer_cprocess.in1process[0][i][last], nex),
			inv_idx(lusifer_cprocess.out1process[0][i][last], nex),
			inv_idx(lusifer_cprocess.virtprocess[0][i][last], nex)
		);

		diagrams.emplace_back(vertices, propagators);
	}

	return diagrams;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template std::vector<diagram> lusifer_diagram_generator(
	std::vector<std::string> const&,
	lusifer_constants<double> const&
);

}
