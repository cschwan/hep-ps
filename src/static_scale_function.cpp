#include "hep/ps/static_scale_function.hpp"

#include <cassert>

namespace hep
{

template <typename T>
static_scale_function<T>::static_scale_function(
	T scale,
	scale_variation variation
)
	: scale_{scale}
	, variation_{variation}
{
}

template <typename T>
void static_scale_function<T>::operator()(
	std::vector<T> const&,
	std::vector<hep::scales<T>>& scales
) {
	switch (variation_)
	{
	case scale_variation::single_scale:
		scales.emplace_back(scale_, scale_, scale_);
		break;

	case scale_variation::three_point:
		scales.emplace_back(         scale_, scale_,          scale_);
		scales.emplace_back(T(0.5) * scale_, scale_, T(0.5) * scale_);
		scales.emplace_back(T(2.0) * scale_, scale_, T(2.0) * scale_);
		break;

	case scale_variation::seven_point:
		scales.emplace_back(         scale_, scale_,          scale_);
		scales.emplace_back(T(0.5) * scale_, scale_, T(0.5) * scale_);
		scales.emplace_back(T(2.0) * scale_, scale_, T(2.0) * scale_);
		scales.emplace_back(T(0.5) * scale_, scale_,          scale_);
		scales.emplace_back(         scale_, scale_, T(0.5) * scale_);
		scales.emplace_back(T(2.0) * scale_, scale_,          scale_);
		scales.emplace_back(         scale_, scale_, T(2.0) * scale_);
		break;

	default:
		assert( false );
	}
}

template <typename T>
bool static_scale_function<T>::dynamic() const
{
	return false;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class static_scale_function<double>;

}
