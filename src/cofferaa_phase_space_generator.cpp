#include "hep/ps/cofferaa_phase_space_generator.hpp"
#include "hep/ps/fortran_helper.hpp"

#include "cofferaa_interfaces.hpp"
#include "hadron_hadron_psg_adapter.hpp"

#include <algorithm>
#include <array>
#include <cassert>

namespace
{

template <typename T>
class cofferaa_psg
{
public:
    using numeric_type = T;

    cofferaa_psg(
        std::vector<int> const& process,
        bool all_diagrams,
        hep::lusifer_constants<T> const& constants,
        std::vector<std::tuple<int, int, int>> const& dipoles
            = std::vector<std::tuple<int, int, int>>{},
        std::size_t extra_random_numbers = 0
    );

    std::size_t channels() const;

    T densities(std::vector<T>& densities);

    std::size_t dimensions() const;

    void generate(
        std::vector<T> const& random_numbers,
        std::vector<T>& momenta,
        T cmf_energy,
        std::size_t channel
    );

    std::size_t map_dimensions() const;

private:
    T cmf_energy_;

    std::size_t channels_;
    std::size_t max_particles_;
    std::size_t particles_;

    std::vector<T> current_point_;
};

template <typename T>
cofferaa_psg<T>::cofferaa_psg(
    std::vector<int> const& process,
    bool all_diagrams,
    hep::lusifer_constants<T> const& constants,
    std::vector<std::tuple<int, int, int>> const& dipoles,
    std::size_t extra_random_numbers
) {
    // TODO: NYI
    assert( extra_random_numbers == 0 );

    particles_ = process.size();

    int maxex;
    int maxgen;
    cofferaa_extra_max(&maxex, &maxgen);

    max_particles_ = maxex;

    assert( max_particles_ >= particles_ );

    current_point_.resize(4 * maxex);

    // FORTRAN counting: use the first generator
    int generator = 1;
    double mw = constants.mass_w;
    double gw = constants.width_w;
    double mz = constants.mass_z;
    double gz = constants.width_z;
    double mh = constants.mass_h;
    double gh = constants.width_h;
    double mt = constants.mass_t;
    double gt = constants.width_t;

    // set constants and the number of particles
    cofferaa_extra_set(
        &generator,
        &mw,
        &gw,
        &mz,
        &gz,
        &mh,
        &gh,
        &mt,
        &gt
    );

    // only used for technical cuts
    double energy = 0.0;
    // output variable; minimum energy that the generator produces
    double smin;
    // number of external particles
    int next = process.size();
    int smodel = all_diagrams ? 12 : 0;
    int sincludecuts = 0;
    // generate phase space according to the dipole mappings
    int ssub = (dipoles.size() != 0) ? 1 : 0;

    std::vector<int> dipole_emitter;
    std::vector<int> dipole_unresolved;
    std::vector<int> dipole_spectator;

    for (auto const& dipole : dipoles)
    {
        dipole_emitter.push_back(std::get<0>(dipole));
        dipole_unresolved.push_back(std::get<1>(dipole));
        dipole_spectator.push_back(std::get<2>(dipole));
    }

    int dipole_count = dipoles.size();

    cofferaa_initgenerator(
        &energy,
        &smin,
        process.data(),
        &generator,
        &next,
        &smodel,
        &sincludecuts,
        &ssub,
        &dipole_count,
        dipole_emitter.data(),
        dipole_unresolved.data(),
        dipole_spectator.data()
    );

    int channels;
    cofferaa_extra_data(&generator, &channels);

    assert( channels > 0 );

    channels_ = channels;
}

template <typename T>
std::size_t cofferaa_psg<T>::channels() const
{
    return channels_;
}

template <typename T>
T cofferaa_psg<T>::densities(std::vector<T>& densities)
{
    using std::acos;
    using std::pow;

    assert( densities.size() >= channels_ );

    int generator = 1;
    int switch_ = 2;

    cofferaa_density(
        current_point_.data(),
        densities.data(),
        &generator,
        &switch_
    );

    return pow(T(0.5) / acos(T(-1.0)), T(3 * particles_ - 10));
}

template <typename T>
std::size_t cofferaa_psg<T>::dimensions() const
{
    return 3 * (particles_ - 4) + 2;
}

template <typename T>
void cofferaa_psg<T>::generate(
    std::vector<T> const& random_numbers,
    std::vector<T>& momenta,
    T cmf_energy,
    std::size_t channel
) {
    assert( momenta.size() == map_dimensions() );

    double const value = 0.5 * static_cast <double> (cmf_energy);
    double kbeam[] = { value, value, 0.0, 0.0, 0.0, 0.0, -value, value };
    // dummy parameter
    double g;
    // FORTRAN counting
    int channel_ = channel + 1;
    int generator = 1;
    int switch_ = 1;

    cofferaa_generation(
        random_numbers.data(),
        &kbeam[0],
        current_point_.data(),
        &g,
        &channel_,
        &generator,
        &switch_
    );

    momenta = current_point_;
    hep::fortran_ordering_to_cpp(momenta);
    momenta.resize(map_dimensions());
}

template <typename T>
std::size_t cofferaa_psg<T>::map_dimensions() const
{
    return 4 * particles_;
}

}

namespace hep
{

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_cofferaa_phase_space_generator(
    T min_energy,
    T cmf_energy,
    std::vector<int> const& process,
    lusifer_constants<T> const& constants,
    std::vector<std::tuple<int, int, int>> const& dipoles,
    std::size_t extra_random_numbers
) {
    return std::make_unique<hadron_hadron_psg_adapter<cofferaa_psg<T>>>(
        min_energy,
        cmf_energy,
        process,
        true,
        constants,
        dipoles,
        extra_random_numbers
    );
}

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_minimal_cofferaa_phase_space_generator(
    T min_energy,
    T cmf_energy,
    std::vector<int> const& process,
    lusifer_constants<T> const& constants,
    std::vector<std::tuple<int, int, int>> const& dipoles,
    std::size_t extra_random_numbers
) {
    return std::make_unique<hadron_hadron_psg_adapter<cofferaa_psg<T>>>(
        min_energy,
        cmf_energy,
        process,
        false,
        constants,
        dipoles,
        extra_random_numbers
    );
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template std::unique_ptr<phase_space_generator<double>>
make_cofferaa_phase_space_generator(
    double,
    double,
    std::vector<int> const&,
    lusifer_constants<double> const&,
    std::vector<std::tuple<int, int, int>> const&,
    std::size_t
);

template std::unique_ptr<phase_space_generator<double>>
make_minimal_cofferaa_phase_space_generator(
    double,
    double,
    std::vector<int> const&,
    lusifer_constants<double> const&,
    std::vector<std::tuple<int, int, int>> const&,
    std::size_t
);

}
