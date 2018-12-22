#include "hep/ps/parton_dfs.hpp"
#include "hep/ps/pdg_functions.hpp"
#include "hep/ps/suppress_banners.hpp"

#include <LHAPDF/PDF.h>

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace hep
{

template <typename T>
class parton_dfs<T>::impl
{
public:
    impl(std::string const& pdf_name, std::size_t pdf_member);
    impl(std::string const& pdf_name);

    std::vector<std::unique_ptr<LHAPDF::PDF>> pdfs;
    std::vector<std::pair<int, parton>> pdg_ids;
};

template <typename T>
parton_dfs<T>::impl::impl(std::string const& pdf_name, std::size_t pdf_member)
{
    if (suppress_banners())
    {
        LHAPDF::setVerbosity(0);
    }

    pdfs.emplace_back(LHAPDF::mkPDF(pdf_name, pdf_member));
}

template <typename T>
parton_dfs<T>::impl::impl(std::string const& pdf_name)
{
    if (suppress_banners())
    {
        LHAPDF::setVerbosity(0);
    }

    for (auto pdf : LHAPDF::mkPDFs(pdf_name))
    {
        pdfs.emplace_back(pdf);
    }
}

template <typename T>
parton_dfs<T>::parton_dfs(
    std::string const& pdf_name,
    std::size_t pdf_member
)
    : pimpl(std::make_unique<impl>(pdf_name, pdf_member))
{
}

template <typename T>
parton_dfs<T>::parton_dfs(std::string const& pdf_name)
    : pimpl(std::make_unique<impl>(pdf_name))
{
}

template <typename T>
parton_dfs<T>::parton_dfs(parton_dfs<T>&& pdf)
    : pimpl(std::move(pdf.pimpl))
{
}

template <typename T>
parton_dfs<T>::~parton_dfs() = default;

template <typename T>
void parton_dfs<T>::eval_alphas(
    nonstd::span<hep::scales<T> const> scales,
    std::vector<T>& alphas
) {
    auto& pdf = pimpl->pdfs.front();

    // TODO: implement a cache?
    for (auto const& scale : scales)
    {
        alphas.push_back(pdf->alphasQ(scale.renormalization()));
    }
}

template <typename T>
std::size_t parton_dfs<T>::count() const
{
    return pimpl->pdfs.size();
}

template <typename T>
void parton_dfs<T>::eval(
    T x,
    std::size_t scale_count,
    nonstd::span<scales<T> const> scales,
    std::vector<parton_array<T>>& scale_pdfs,
    std::vector<parton_array<T>>& uncertainty_pdfs
) {
    scale_pdfs.clear();
    scale_pdfs.resize(scales.size());

    auto const xval = static_cast <double> (x);
    std::size_t const pdf_count = count();
    std::size_t const central_scales = scales.size() / scale_count;

    if (pdf_count > 1)
    {
        uncertainty_pdfs.clear();
        uncertainty_pdfs.resize(central_scales * pdf_count);
    }

    for (std::size_t i = 0; i != scales.size(); ++i)
    {
        // check if we already calculated the central PDF for the current factorization scale
        auto const end = std::next(scales.begin(), i);
        auto const result = std::find_if(scales.begin(), end, [&](hep::scales<T> const& s) {
            return s.factorization() == scales[i].factorization();
        });
        std::size_t const index = std::distance(scales.begin(), result);

        auto const muf = static_cast <double> (scales[i].factorization());

        if (result == end)
        {
            for (auto const& ids : pimpl->pdg_ids)
            {
                auto const xfx = pimpl->pdfs.front()->xfxQ(ids.first, xval, muf);
                scale_pdfs.at(i)[ids.second] = T(xfx) / x;
            }
        }
        else
        {
            scale_pdfs.at(i) = scale_pdfs.at(index);
        }

        // check if we are evaluating PDFs for the central scale
        if ((pdf_count > 1) && (i % central_scales == 0))
        {
            std::size_t const j = i / scale_count;

            // we already calculated the central PDFs
            uncertainty_pdfs.at(j * pdf_count) = scale_pdfs.at(i);

            if (result == end)
            {
                for (std::size_t k = 1; k != pdf_count; ++k)
                {
                    for (auto const& ids : pimpl->pdg_ids)
                    {
                        auto const xfx = pimpl->pdfs.front()->xfxQ(ids.first, xval, muf);
                        uncertainty_pdfs.at(j * pdf_count + k)[ids.second] = T(xfx) / x;
                    }
                }
            }
            else
            {
                uncertainty_pdfs.insert(
                    uncertainty_pdfs.end(),
                    std::next(uncertainty_pdfs.begin(), pdf_count * index),
                    std::next(uncertainty_pdfs.begin(), pdf_count * (index + 1))
                );
            }
        }
    }
}

template <typename T>
void parton_dfs<T>::register_partons(parton_set set)
{
    for (auto const p : set)
    {
        int const pdg_id = parton_to_pdg_id(p);

        if (!pimpl->pdfs.front()->hasFlavor(pdg_id))
        {
            throw std::domain_error("PDF doesn't support the requested flavor");
        }

        pimpl->pdg_ids.emplace_back(pdg_id, p);
    }

    pimpl->pdg_ids.shrink_to_fit();
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class parton_dfs<double>;

}
