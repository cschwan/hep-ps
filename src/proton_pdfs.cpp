#include "hep/ps/proton_pdfs.hpp"

#include <LHAPDF/PDF.h>

#include <utility>

namespace hep
{

template <typename T>
class proton_pdfs<T>::impl
{
public:
	impl(std::string const& pdf_name, std::size_t pdf_member);
	impl(std::string const& pdf_name);

	std::vector<std::unique_ptr<LHAPDF::PDF>> pdfs;
	std::vector<double> xfx;
	alphas_calc<T> alphas;
};

template <typename T>
proton_pdfs<T>::impl::impl(std::string const& pdf_name, std::size_t pdf_member)
	: pdfs()
	// TODO: is it possible to extract this number from LHAPDF?
	, xfx(13)
	, alphas(nullptr)
{
	pdfs.emplace_back(LHAPDF::mkPDF(pdf_name, pdf_member));
	alphas = alphas_calc<T>(pdfs.front().get());
}

// TODO: fix the implementation
template <typename T>
proton_pdfs<T>::impl::impl(std::string const& pdf_name)
	: pdfs()
	// TODO: is it possible to extract this number from LHAPDF?
	, xfx(13)
	, alphas(nullptr)
{
	for (auto pdf : LHAPDF::mkPDFs(pdf_name))
	{
		pdfs.emplace_back(pdf);
	}

	alphas = alphas_calc<T>(pdfs.front().get());
}

template <typename T>
proton_pdfs<T>::proton_pdfs(
	std::string const& pdf_name,
	std::size_t pdf_member
)
	// TODO: replace with `make_unique` in C++14
	: pimpl(new impl(pdf_name, pdf_member))
{
}

template <typename T>
proton_pdfs<T>::proton_pdfs(std::string const& pdf_name)
	: pimpl(new impl(pdf_name))
{
}

template <typename T>
proton_pdfs<T>::proton_pdfs(proton_pdfs<T>&& pdf)
	: pimpl(std::move(pdf.pimpl))
{
}

template <typename T>
proton_pdfs<T>::~proton_pdfs() = default;

template <typename T>
alphas_calc<T>& proton_pdfs<T>::alphas()
{
	return pimpl->alphas;
}

template <typename T>
std::size_t proton_pdfs<T>::count() const
{
	return pimpl->pdfs.size();
}

template <typename T>
void proton_pdfs<T>::eval(T x, T scale, std::vector<parton_array<T>>& pdfs)
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

template class proton_pdfs<double>;

}
