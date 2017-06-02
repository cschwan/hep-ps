#include "hep/ps/lusifer_diagram_generator.hpp"

#include "lusifer_interfaces.hpp"

#include <algorithm>
#include <cassert>

constexpr std::size_t maxg = 1;
constexpr std::size_t maxe = 9;
constexpr std::size_t maxch = 20000;
constexpr std::size_t maxv = 40;

extern "C"
{

extern struct
{
	double alphaisr;        // NOT NEEDED
	double scale;           // NOT NEEDED
	double meisr;           // NOT NEEDED
	double s[1 << maxe];    // NOT NEEDED
	double p[1 << maxe][4]; // NOT NEEDED
	double mass[maxv + 1];
	double width[maxv + 1];
	int nchannel[maxg];
	int nexternal[maxg];    // NOT NEEDED
	int allbinary[maxg];    // NOT NEEDED
} lusifer_general_;

extern struct
{
	double powerinv[maxg][maxch][maxe];    // NOT NEEDED
	double mcutinv[maxg][(1 << maxe) + 1];
	int ininv[maxg][maxch][maxe];
	int idhepinv[maxg][maxch][maxv];
	int ninv[maxg][maxch];
	int lmin[maxg][maxch][maxe][maxe];
	int lmax[maxg][maxch][maxe][maxe];
} lusifer_cinv_;

extern struct
{
	int indecay[maxg][maxch][maxe];
	int out1decay[maxg][maxch][maxe];
	int out2decay[maxg][maxch][maxe];
	int ndecay[maxg][maxch];
} lusifer_cdecay_;

extern struct
{
	double powerprocess[maxg][maxch][maxe];
	double ccutprocess[maxg][1 << maxe];
	int in1process[maxg][maxch][maxe];
	int in2process[maxg][maxch][maxe];
	int out1process[maxg][maxch][maxe];
	int out2process[maxg][maxch][maxe];
	int inprocess[maxg][maxch][maxe];
	int virtprocess[maxg][maxch][maxe];
	int idhepprocess[maxg][maxch][maxe];
	int nprocess[maxg][maxch];
} lusifer_cprocess_;

extern struct
{
	int nsinv[maxg][maxch];
	int chinv[maxg][maxch];
	int maxinv[maxg];
	int ntprocess[maxg][maxch];
	int chprocess[maxg][maxch];
	int maxprocess[maxg];
	int nsdecay[maxg][maxch];
	int chdecay[maxg][maxch];
	int maxdecay[maxg];
	int numinv[maxg][maxch][maxe];
	int numprocess[maxg][maxch][maxe];
	int numdecay[maxg][maxch][maxe];
} lusifer_cdensity_;

}

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
