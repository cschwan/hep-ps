#include "hep/ps/lusifer_diagram_generator.hpp"

#include "lusifer_interfaces.hpp"

#include <algorithm>
#include <cassert>

namespace hep
{

template <typename T>
std::vector<diagram> lusifer_diagram_generator(
	std::vector<std::string> const& processes,
	lusifer_constants<T> const& constants
) {
	int maxex;
	int maxgen;
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
			lusifer_extra_set(&g, &nex, &mw, &gw, &mz, &gz, &mh, &gh, &mt, &gt);

			// the number of particles must be supported by the generator
			assert( nex <= maxex );
		}

		// all processes must have the same number of particles
		assert( nex == static_cast <int> (process.size() / 3) );

		// fill up the string with three spaces for particle not used
		std::string process0 = process;
		process0.append(3 * (maxex - nex), ' ');

		// TODO: what is the meaning of this parameter?
		int lightfermions = 0;
		// do not include cuts in the phase space generation
		int includecuts = 0;
		// do not print channel information
		int sout = 0;

		lusifer_initphasespace(
			process0.c_str(),
			&g,
			&lightfermions,
			&includecuts,
			&sout,
			process0.size()
		);
	}

	int channels;
	lusifer_extra_data(&g, &channels);

	// there must be at least one channel, otherwise something went wrong
	assert( channels > 0 );

	std::vector<diagram> diagrams;
	diagrams.reserve(channels);

	for (int i = 0; i != channels; ++i)
	{
		std::vector<vertex> vertices;
		std::vector<propagator> propagators;

		for (int j = 0; j != lusifer_cdecay_.ndecay[0][i]; ++j)
		{
			vertices.emplace_back(
				lusifer_cdecay_.indecay[0][i][j],
				lusifer_cdecay_.out1decay[0][i][j],
				lusifer_cdecay_.out2decay[0][i][j]
			);

			int k = 0;
			for (; k != lusifer_cinv_.ninv[0][i]; ++k)
			{
				if (lusifer_cdecay_.indecay[0][i][j] ==
					lusifer_cinv_.ininv[0][i][k])
				{
					break;
				}
			}

			assert( k != lusifer_cinv_.ninv[0][i] );

			propagators.emplace_back(
				lusifer_cdecay_.indecay[0][i][j],
				lusifer_cinv_.idhepinv[0][i][k]
			);
		}

		for (int j = 0; j != lusifer_cprocess_.nprocess[0][i]; ++j)
		{
			vertices.emplace_back(
				lusifer_cprocess_.in1process[0][i][j],
				lusifer_cprocess_.out1process[0][i][j],
				lusifer_cprocess_.virtprocess[0][i][j]
			);

			propagators.emplace_back(
				lusifer_cprocess_.idhepprocess[0][i][j],
				lusifer_cprocess_.virtprocess[0][i][j]
			);
		}

		int last = lusifer_cprocess_.nprocess[0][i] - 1;

		vertices.emplace_back(
			lusifer_cprocess_.in1process[0][i][last],
			lusifer_cprocess_.out1process[0][i][last],
			lusifer_cprocess_.virtprocess[0][i][last]
		);

		if (lusifer_cprocess_.nprocess[0][i] == 0)
		{
			// TODO: if there are no t-channel propagators, add initial state
			// s-channel propagator
		}

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
