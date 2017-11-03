#include "hep/ps/parton_dfs.hpp"

#include <LHAPDF/PDF.h>

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
	std::vector<hep::scales<T>> const& scales,
	std::vector<parton_array<T>>& pdfs
) {
	// TODO: is it possible to get these values from LHAPDF?
	constexpr std::size_t lhapdf_ac = -4 + 6;
	constexpr std::size_t lhapdf_as = -3 + 6;
	constexpr std::size_t lhapdf_au = -2 + 6;
	constexpr std::size_t lhapdf_ad = -1 + 6;
	constexpr std::size_t lhapdf_g  =  0 + 6;
	constexpr std::size_t lhapdf_d  =  1 + 6;
	constexpr std::size_t lhapdf_u  =  2 + 6;
	constexpr std::size_t lhapdf_s  =  3 + 6;
	constexpr std::size_t lhapdf_c  =  4 + 6;

	// TODO: implement a cache?
	for (std::size_t i = 0; i != scales.size(); ++i)
	{
		pimpl->pdfs.front()->xfxQ(
			static_cast <double> (x),
			static_cast <double> (scales.at(i).factorization()),
			pimpl->xfx
		);

		pdfs.at(i)[parton::anti_charm]   = T(pimpl->xfx[lhapdf_ac]) / x;
		pdfs.at(i)[parton::anti_strange] = T(pimpl->xfx[lhapdf_as]) / x;
		pdfs.at(i)[parton::anti_up]      = T(pimpl->xfx[lhapdf_au]) / x;
		pdfs.at(i)[parton::anti_down]    = T(pimpl->xfx[lhapdf_ad]) / x;
		pdfs.at(i)[parton::gluon]        = T(pimpl->xfx[lhapdf_g])  / x;
		pdfs.at(i)[parton::down]         = T(pimpl->xfx[lhapdf_d])  / x;
		pdfs.at(i)[parton::up]           = T(pimpl->xfx[lhapdf_u])  / x;
		pdfs.at(i)[parton::strange]      = T(pimpl->xfx[lhapdf_s])  / x;
		pdfs.at(i)[parton::charm]        = T(pimpl->xfx[lhapdf_c])  / x;
	}
}

template <typename T>
void parton_dfs<T>::eval(T x, T scale, std::vector<parton_array<T>>& pdfs)
{
	// TODO: is it possible to get these values from LHAPDF?
	constexpr std::size_t lhapdf_ac = -4 + 6;
	constexpr std::size_t lhapdf_as = -3 + 6;
	constexpr std::size_t lhapdf_au = -2 + 6;
	constexpr std::size_t lhapdf_ad = -1 + 6;
	constexpr std::size_t lhapdf_g  =  0 + 6;
	constexpr std::size_t lhapdf_d  =  1 + 6;
	constexpr std::size_t lhapdf_u  =  2 + 6;
	constexpr std::size_t lhapdf_s  =  3 + 6;
	constexpr std::size_t lhapdf_c  =  4 + 6;

	for (std::size_t i = 0; i != pdfs.size(); ++i)
	{
		pimpl->pdfs.at(i)->xfxQ(
			static_cast <double> (x),
			static_cast <double> (scale),
			pimpl->xfx
		);

		pdfs.at(i)[parton::anti_charm]   = T(pimpl->xfx[lhapdf_ac]) / x;
		pdfs.at(i)[parton::anti_strange] = T(pimpl->xfx[lhapdf_as]) / x;
		pdfs.at(i)[parton::anti_up]      = T(pimpl->xfx[lhapdf_au]) / x;
		pdfs.at(i)[parton::anti_down]    = T(pimpl->xfx[lhapdf_ad]) / x;
		pdfs.at(i)[parton::gluon]        = T(pimpl->xfx[lhapdf_g])  / x;
		pdfs.at(i)[parton::down]         = T(pimpl->xfx[lhapdf_d])  / x;
		pdfs.at(i)[parton::up]           = T(pimpl->xfx[lhapdf_u])  / x;
		pdfs.at(i)[parton::strange]      = T(pimpl->xfx[lhapdf_s])  / x;
		pdfs.at(i)[parton::charm]        = T(pimpl->xfx[lhapdf_c])  / x;
	}
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class parton_dfs<double>;

}
