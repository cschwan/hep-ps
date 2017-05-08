#include "hep/ps/lusifer_constants.hpp"

namespace hep
{

template <typename T>
lusifer_constants<T>::lusifer_constants(
	T mass_h,
	T width_h,
	T mass_t,
	T width_t,
	T mass_w,
	T width_w,
	T mass_z,
	T width_z
)
	: mass_h(mass_h)
	, width_h(width_h)
	, mass_t(mass_t)
	, width_t(width_t)
	, mass_w(mass_w)
	, width_w(width_w)
	, mass_z(mass_z)
	, width_z(width_z)
{
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class lusifer_constants<double>;

}
