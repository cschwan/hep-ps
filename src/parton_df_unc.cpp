#include "hep/ps/parton_df_unc.hpp"

#include <LHAPDF/PDFSet.h>

#include <cmath>
#include <utility>

namespace hep
{

template <typename T>
class parton_df_unc<T>::impl
{
public:
	impl(std::string const& name, T cl);

	LHAPDF::PDFSet pdfset;
	T cl;
};

template <typename T>
parton_df_unc<T>::impl::impl(std::string const& name, T cl)
	: pdfset{name}
	, cl{cl}
{
	using std::erf;
	using std::sqrt;

	if (cl == T())
	{
		this->cl = T(100.0) * erf(T(1.0) / sqrt(T(2.0)));
	}
}

template <typename T>
parton_df_unc<T>::parton_df_unc(std::string const& name, T cl)
	: pimpl(std::make_unique<impl>(name, cl))
{
}

template <typename T>
parton_df_unc<T>::parton_df_unc(parton_df_unc<T>&& pdf_unc)
	: pimpl(std::move(pdf_unc.pimpl))
{
}

template <typename T>
parton_df_unc<T>::~parton_df_unc() = default;

template <typename T>
pdf_unc<T> parton_df_unc<T>::uncertainty(std::vector<T> const& values) const
{
	std::vector<double> double_values;
	double_values.reserve(values.size());

	for (T const value : values)
	{
		double_values.push_back(static_cast <double> (value));
	}

	auto const unc = pimpl->pdfset.uncertainty(
		double_values,
		static_cast <double> (pimpl->cl)
	);

	pdf_unc<T> result;
	result.sym = T(unc.errsymm);
	result.neg = T(unc.errminus);
	result.pos = T(unc.errplus);

	return result;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class parton_df_unc<double>;

}
