#include "hep/ps/proton_pdf.hpp"

#include <LHAPDF/PDF.h>

#include <utility>
#include <vector>

namespace hep
{

template <typename T>
class proton_pdf<T>::impl
{
public:
	impl(std::string const& pdf_name, std::size_t pdf_member);

	std::unique_ptr<LHAPDF::PDF> pdf;
	std::vector<double> xfx;
	alphas_calc<T> alphas;
};

template <typename T>
proton_pdf<T>::impl::impl(
	std::string const& pdf_name,
	std::size_t pdf_member
)
	: pdf(LHAPDF::mkPDF(pdf_name, pdf_member))
	// TODO: is it possible to extract this number from LHAPDF?
	, xfx(13)
	, alphas(pdf.get())
{
}

template <typename T>
proton_pdf<T>::proton_pdf(
	std::string const& pdf_name,
	std::size_t pdf_member
)
	// TODO: replace with `make_unique` in C++14
	: pimpl(new impl(pdf_name, pdf_member))
{
}

template <typename T>
proton_pdf<T>::proton_pdf(proton_pdf<T>&& pdf)
	: pimpl(std::move(pdf.pimpl))
{
}

template <typename T>
proton_pdf<T>::~proton_pdf() = default;

template <typename T>
alphas_calc<T>& proton_pdf<T>::alphas()
{
	return pimpl->alphas;
}

template <typename T>
parton_array<T> proton_pdf<T>::pdf(T x, T scale)
{
	parton_array<T> result;

	pimpl->pdf->xfxQ(
		static_cast <double> (x),
		static_cast <double> (scale),
		pimpl->xfx
	);

	// TODO: is it possible to get these values from LHAPDF?
	constexpr std::size_t lhapdf_anti_charm   = -4 + 6;
	constexpr std::size_t lhapdf_anti_strange = -3 + 6;
	constexpr std::size_t lhapdf_anti_up      = -2 + 6;
	constexpr std::size_t lhapdf_anti_down    = -1 + 6;
	constexpr std::size_t lhapdf_gluon        =  0 + 6;
	constexpr std::size_t lhapdf_down         =  1 + 6;
	constexpr std::size_t lhapdf_up           =  2 + 6;
	constexpr std::size_t lhapdf_strange      =  3 + 6;
	constexpr std::size_t lhapdf_charm        =  4 + 6;

	result[parton::anti_charm]   = T(pimpl->xfx[lhapdf_anti_charm])   / x;
	result[parton::anti_strange] = T(pimpl->xfx[lhapdf_anti_strange]) / x;
	result[parton::anti_up]      = T(pimpl->xfx[lhapdf_anti_up])      / x;
	result[parton::anti_down]    = T(pimpl->xfx[lhapdf_anti_down])    / x;
	result[parton::gluon]        = T(pimpl->xfx[lhapdf_gluon])        / x;
	result[parton::down]         = T(pimpl->xfx[lhapdf_down])         / x;
	result[parton::up]           = T(pimpl->xfx[lhapdf_up])           / x;
	result[parton::strange]      = T(pimpl->xfx[lhapdf_strange])      / x;
	result[parton::charm]        = T(pimpl->xfx[lhapdf_charm])        / x;

	return result;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class proton_pdf<double>;

}
