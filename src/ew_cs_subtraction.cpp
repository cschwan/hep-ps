#include "hep/ps/ew_cs_subtraction.hpp"

#include <cmath>

namespace hep
{

using std::sqrt;

template <typename T>
ew_cs_subtraction<T>::ew_cs_subtraction(
	T nf,
	factorization_scheme fscheme,
	renormalization_scheme rscheme
)
	: cs_subtraction<T>(
		T(0.5) * (T(1.0) + sqrt(T(5.0))), // this value of nc sets Cf to one
		T(1.0),
		nf,
		fscheme,
		rscheme
	  )
{
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class ew_cs_subtraction<double>;

}
