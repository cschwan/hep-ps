#include "hep/ps/lusifer_ps_channels.hpp"

#include "lusifer_interfaces.hpp"

#include <cassert>
#include <utility>

namespace hep
{

// TODO: remove dependence on a numerical type `T`

template <typename T>
std::vector<ps_channel> lusifer_ps_channels(
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

            double mw = static_cast <double> (constants.mass_w);
            double gw = static_cast <double> (constants.width_w);
            double mz = static_cast <double> (constants.mass_z);
            double gz = static_cast <double> (constants.width_z);
            double mh = static_cast <double> (constants.mass_h);
            double gh = static_cast <double> (constants.width_h);
            double mt = static_cast <double> (constants.mass_t);
            double gt = static_cast <double> (constants.width_t);

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

        // TODO: what is the meaning of this parameter?
        int lightfermions = 0;
        // do not include cuts in the phase space generation
        int includecuts = 0;

        lusifer_initphasespace(process0.c_str(), g, lightfermions, includecuts);
    }

    int channel_count;
    lusifer_extra_data(g, &channel_count);

    // there must be at least one channel, otherwise something went wrong
    assert( channel_count > 0 );

    std::vector<ps_channel> channels;
    channels.reserve(channel_count);

    for (std::size_t i = 0; i != channels.capacity(); ++i)
    {
        std::vector<ps_invariant> invariants;
        std::vector<ps_tchannel> tchannels;
        std::vector<ps_decay> decays;

        invariants.reserve(lusifer_cinv.ninv[0][i]);
        tchannels.reserve(lusifer_cprocess.nprocess[0][i]);
        decays.reserve(lusifer_cdecay.ndecay[0][i]);

        for (std::size_t j = 0; j != invariants.capacity(); ++j)
        {
            invariants.emplace_back(
                lusifer_cinv.ininv[0][i][j] - 1,
                lusifer_cinv.idhepinv[0][i][j]
            );
        }

        for (std::size_t j = 0; j != tchannels.capacity(); ++j)
        {
            tchannels.emplace_back(
                lusifer_cprocess.in1process[0][i][j] - 1,
                lusifer_cprocess.in2process[0][i][j] - 1,
                lusifer_cprocess.out1process[0][i][j] - 1,
                lusifer_cprocess.out2process[0][i][j] - 1,
                lusifer_cprocess.inprocess[0][i][j] - 1,
                lusifer_cprocess.virtprocess[0][i][j] - 1,
                lusifer_cprocess.idhepprocess[0][i][j]
            );
        }

        for (std::size_t j = 0; j != decays.capacity(); ++j)
        {
            decays.emplace_back(
                lusifer_cdecay.indecay[0][i][j] - 1,
                lusifer_cdecay.out1decay[0][i][j] - 1,
                lusifer_cdecay.out2decay[0][i][j] - 1
            );
        }

        channels.emplace_back(std::move(invariants), std::move(tchannels), std::move(decays));
    }

    return channels;
}

template <typename T>
std::vector<ps_channel> lusifer_ps_channels(
    std::string const& process,
    lusifer_constants<T> const& constants
) {
    return lusifer_ps_channels(std::vector<std::string>{process}, constants);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template std::vector<ps_channel> lusifer_ps_channels(
    std::vector<std::string> const&,
    lusifer_constants<double> const&
);

template std::vector<ps_channel> lusifer_ps_channels(
    std::string const&,
    lusifer_constants<double> const&
);

}
