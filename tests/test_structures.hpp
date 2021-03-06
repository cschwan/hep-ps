#ifndef TEST_STRUCTURES_HPP
#define TEST_STRUCTURES_HPP

#include "hep/mc/distribution_parameters.hpp"

#include "hep/ps/final_state.hpp"
#include "hep/ps/initial_state.hpp"
#include "hep/ps/parton.hpp"
#include "hep/ps/psp.hpp"
#include "hep/ps/recombined_state.hpp"
#include "hep/ps/scales.hpp"

#include "catch2/catch.hpp"

#include "nonstd/span.hpp"

#include <array>
#include <cstddef>
#include <numeric>
#include <vector>

namespace
{

template <typename T>
class test_pdf
{
public:
    void eval(
        T x,
        std::vector<hep::scales<T>> const& scales,
        std::vector<hep::parton_array<T>>& pdfs
    ) {
        CHECK( x >= T() );
        CHECK( x < T(1.0) );
        CHECK( scales.size() == 1 );
        CHECK( pdfs.size() >= 1 );
        CHECK( scales.front().factorization() > T() );

        pdfs.front()[hep::parton::anti_up]      = T(1.0);
        pdfs.front()[hep::parton::anti_down]    = T(1.0);
        pdfs.front()[hep::parton::anti_charm]   = T(1.0);
        pdfs.front()[hep::parton::anti_strange] = T(1.0);
        pdfs.front()[hep::parton::gluon]        = T(1.0);
        pdfs.front()[hep::parton::up]           = T(1.0);
        pdfs.front()[hep::parton::down]         = T(1.0);
        pdfs.front()[hep::parton::charm]        = T(1.0);
        pdfs.front()[hep::parton::strange]      = T(1.0);
    }

    void eval(T x, T scale, std::vector<hep::parton_array<T>>& pdfs)
    {
        CHECK( pdfs.size() == count() );
        CHECK( x >= T() );
        CHECK( x < T(1.0) );
        CHECK( scale > T() );

        pdfs.front()[hep::parton::anti_up]      = T(1.0);
        pdfs.front()[hep::parton::anti_down]    = T(1.0);
        pdfs.front()[hep::parton::anti_charm]   = T(1.0);
        pdfs.front()[hep::parton::anti_strange] = T(1.0);
        pdfs.front()[hep::parton::gluon]        = T(1.0);
        pdfs.front()[hep::parton::up]           = T(1.0);
        pdfs.front()[hep::parton::down]         = T(1.0);
        pdfs.front()[hep::parton::charm]        = T(1.0);
        pdfs.front()[hep::parton::strange]      = T(1.0);
    }

    void eval(
        T x,
        std::size_t /*scale_count*/,
        nonstd::span<hep::scales<T> const> scales,
        std::vector<hep::parton_array<T>>& scale_pdfs,
        std::vector<hep::parton_array<T>>& uncertainty_pdfs
    ) {
        CHECK( x >= T() );
        CHECK( x < T(1.0) );
        CHECK( scales.size() >= 1 );
        CHECK( scales[0].factorization() > T() );
        CHECK( uncertainty_pdfs.empty() );

        scale_pdfs.clear();
        scale_pdfs.resize(scales.size());

        for (std::size_t i = 0; i != size(scales); ++i)
        {
            scale_pdfs.at(i)[hep::parton::anti_up]      = T(1.0);
            scale_pdfs.at(i)[hep::parton::anti_down]    = T(1.0);
            scale_pdfs.at(i)[hep::parton::anti_charm]   = T(1.0);
            scale_pdfs.at(i)[hep::parton::anti_strange] = T(1.0);
            scale_pdfs.at(i)[hep::parton::gluon]        = T(1.0);
            scale_pdfs.at(i)[hep::parton::up]           = T(1.0);
            scale_pdfs.at(i)[hep::parton::down]         = T(1.0);
            scale_pdfs.at(i)[hep::parton::charm]        = T(1.0);
            scale_pdfs.at(i)[hep::parton::strange]      = T(1.0);
        }
    }

    T eval_alphas(T scale)
    {
        CHECK( scale > T() );

        return T(1.0);
    }

    void eval_alphas(nonstd::span<hep::scales<T> const> scales, std::vector<T>& alphas)
    {
        CHECK( scales.size() >= 1 );
        CHECK( scales[0].renormalization() > T() );

        alphas.assign(scales.size(), T(1.0));
    }

    std::size_t count() const
    {
        return 1;
    }

    void register_partons(hep::parton_set)
    {
    }
};

template <typename T>
class test_matrix_elements
{
public:
    test_matrix_elements(hep::initial_state_set set, std::size_t final_states)
        : set_(set)
        , final_states_(final_states)
    {
    }

    void borns(
        std::vector<T> const& phase_space,
        hep::initial_state_set set,
        nonstd::span<hep::scales<T> const> scales,
        std::vector<hep::initial_state_map<T>>& me
    ) const {
        CHECK( phase_space.size() == 4 * (final_states_ + 2) );
        CHECK( set == set_ );
        CHECK( size(scales) == me.size() );

        for (auto const process : set)
        {
            for (std::size_t i = 0; i != size(scales); ++i)
            {
                me.at(i).emplace_back(process, T(1.0));
            }
        }
    }

    std::vector<hep::final_state> final_states() const
    {
        return std::vector<hep::final_state>(final_states_,
            hep::final_state::quark_gluon);
    }

    void alphas(T)
    {
    }

    std::size_t alphas_power() const
    {
        return 0;
    }

private:
    hep::initial_state_set set_;
    std::size_t final_states_;
};

template <typename T>
class test_scale_setter
{
public:
    test_scale_setter(bool dynamic)
        : scale{T(1.0)}
        , dynamic_{dynamic}
    {
    }

    void eval(hep::psp<T> const&, nonstd::span<hep::scales<T>> scales)
    {
        scales[0] = hep::scales<T>{scale, scale, scale};
    }

    bool dynamic() const
    {
        return dynamic_;
    }

    std::size_t count() const
    {
        return 1;
    }

private:
    T scale;
    bool dynamic_;
};

template <typename T>
inline std::vector<hep::distribution_parameters<T>>
test_distribution_parameters()
{
    return std::vector<hep::distribution_parameters<T>>();
}

}

#endif
