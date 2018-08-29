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
    std::vector<hep::scales<T>> const& scales,
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
    std::vector<scales<T>> const& scales,
    std::vector<parton_array<T>>& scale_pdfs,
    std::vector<parton_array<T>>& uncertainty_pdfs
) {
    scale_pdfs.clear();
    scale_pdfs.resize(scales.size());

    auto const xval = static_cast <double> (x);

    for (std::size_t i = 0; i != scales.size(); ++i)
    {
        // check if we already calculated the PDF
        auto const end = std::next(scales.begin(), i);
        auto const result = std::find_if(scales.begin(), end,
            [&](hep::scales<T> const& s) {
                return s.factorization() == scales.at(i).factorization();
        });

        if (result == end)
        {
            auto const muf =
                static_cast <double> (scales.at(i).factorization());

            for (auto const& ids : pimpl->pdg_ids)
            {
                auto const xfx = pimpl->pdfs.front()->xfxQ(
                    ids.first,
                    xval,
                    muf
                );

                scale_pdfs.at(i)[ids.second] = T(xfx) / x;
            }
        }
        else
        {
            std::size_t const index = std::distance(scales.begin(), result);
            scale_pdfs.at(i) = scale_pdfs.at(index);
        }
    }

    std::size_t const size = count();

    if (size > 1)
    {
        uncertainty_pdfs.clear();
        uncertainty_pdfs.resize(size);

        // we already calculated the central PDF
        uncertainty_pdfs.front() = scale_pdfs.front();
    }

    auto const muf = static_cast <double> (scales.front().factorization());

    for (std::size_t i = 1; i != size; ++i)
    {
        for (auto const& ids : pimpl->pdg_ids)
        {
            auto const xfx = pimpl->pdfs.front()->xfxQ(ids.first, xval, muf);
            uncertainty_pdfs.at(i)[ids.second] = T(xfx) / x;
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
