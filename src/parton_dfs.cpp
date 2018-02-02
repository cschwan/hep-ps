#include "hep/ps/parton_dfs.hpp"

#include <LHAPDF/PDF.h>

#include <algorithm>
#include <iterator>
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
	std::vector<double> xfx;
};

template <typename T>
parton_dfs<T>::impl::impl(std::string const& pdf_name, std::size_t pdf_member)
	: pdfs()
	// TODO: is it possible to extract this number from LHAPDF?
	, xfx(13)
{
	pdfs.emplace_back(LHAPDF::mkPDF(pdf_name, pdf_member));
}

// TODO: fix the implementation
template <typename T>
parton_dfs<T>::impl::impl(std::string const& pdf_name)
	: pdfs()
	// TODO: is it possible to extract this number from LHAPDF?
	, xfx(13)
{
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
	// TODO: is it possible to get these values from LHAPDF?
	constexpr std::size_t cx = -4 + 6;
	constexpr std::size_t sx = -3 + 6;
	constexpr std::size_t ux = -2 + 6;
	constexpr std::size_t dx = -1 + 6;
	constexpr std::size_t gl =  0 + 6;
	constexpr std::size_t dq =  1 + 6;
	constexpr std::size_t uq =  2 + 6;
	constexpr std::size_t sq =  3 + 6;
	constexpr std::size_t cq =  4 + 6;

	scale_pdfs.clear();
	scale_pdfs.resize(scales.size());

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
			pimpl->pdfs.front()->xfxQ(
				static_cast <double> (x),
				static_cast <double> (scales.at(i).factorization()),
				pimpl->xfx
			);

			scale_pdfs.at(i)[parton::anti_charm]   = T(pimpl->xfx[cx]) / x;
			scale_pdfs.at(i)[parton::anti_strange] = T(pimpl->xfx[sx]) / x;
			scale_pdfs.at(i)[parton::anti_up]      = T(pimpl->xfx[ux]) / x;
			scale_pdfs.at(i)[parton::anti_down]    = T(pimpl->xfx[dx]) / x;
			scale_pdfs.at(i)[parton::gluon]        = T(pimpl->xfx[gl]) / x;
			scale_pdfs.at(i)[parton::down]         = T(pimpl->xfx[dq]) / x;
			scale_pdfs.at(i)[parton::up]           = T(pimpl->xfx[uq]) / x;
			scale_pdfs.at(i)[parton::strange]      = T(pimpl->xfx[sq]) / x;
			scale_pdfs.at(i)[parton::charm]        = T(pimpl->xfx[cq]) / x;
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

	for (std::size_t i = 1; i != size; ++i)
	{
		pimpl->pdfs.front()->xfxQ(
			static_cast <double> (x),
			static_cast <double> (scales.front().factorization()),
			pimpl->xfx
		);

		uncertainty_pdfs.at(i)[parton::anti_charm]   = T(pimpl->xfx[cx]) / x;
		uncertainty_pdfs.at(i)[parton::anti_strange] = T(pimpl->xfx[sx]) / x;
		uncertainty_pdfs.at(i)[parton::anti_up]      = T(pimpl->xfx[ux]) / x;
		uncertainty_pdfs.at(i)[parton::anti_down]    = T(pimpl->xfx[dx]) / x;
		uncertainty_pdfs.at(i)[parton::gluon]        = T(pimpl->xfx[gl]) / x;
		uncertainty_pdfs.at(i)[parton::down]         = T(pimpl->xfx[dq]) / x;
		uncertainty_pdfs.at(i)[parton::up]           = T(pimpl->xfx[uq]) / x;
		uncertainty_pdfs.at(i)[parton::strange]      = T(pimpl->xfx[sq]) / x;
		uncertainty_pdfs.at(i)[parton::charm]        = T(pimpl->xfx[cq]) / x;
	}
}

template <typename T>
void parton_dfs<T>::register_partons(parton_set set)
{
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class parton_dfs<double>;

}
