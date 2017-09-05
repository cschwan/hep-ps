#include "hep/ps/alphas_calc.hpp"

#include <cassert>

#include <LHAPDF/PDF.h>

namespace hep
{

template <typename T>
alphas_calc<T>::alphas_calc(void* pdf)
	: pdf_(pdf)
{
}

template <typename T>
alphas_calc<T>& alphas_calc<T>::operator=(alphas_calc<T>&& other)
{
	pdf_ = other.pdf_;
	other.pdf_ = nullptr;

	return *this;
}

template <typename T>
T alphas_calc<T>::alphas(T scale)
{
	assert( pdf_ != nullptr );

	return static_cast <LHAPDF::PDF*> (pdf_)->alphasQ(scale);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class alphas_calc<double>;

}
