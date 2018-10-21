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
    std::vector<hep::scales<T>>& scales,
    std::vector<recombined_state> const&
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
void static_scale_function<T>::operator()(psp<T> const&, std::vector<hep::scales<T>>& scales)
{
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
void static_scale_function<T>::eval(psp<T> const&, nonstd::span<hep::scales<T>> scales)
{
    switch (variation_)
    {
    case scale_variation::single_scale:
        scales[0] = hep::scales<T>{scale_, scale_, scale_};
        break;

    case scale_variation::three_point:
        scales[0] = hep::scales<T>{         scale_, scale_,          scale_};
        scales[1] = hep::scales<T>{T(0.5) * scale_, scale_, T(0.5) * scale_};
        scales[2] = hep::scales<T>{T(2.0) * scale_, scale_, T(2.0) * scale_};
        break;

    case scale_variation::seven_point:
        scales[0] = hep::scales<T>{         scale_, scale_,          scale_};
        scales[0] = hep::scales<T>{T(0.5) * scale_, scale_, T(0.5) * scale_};
        scales[0] = hep::scales<T>{T(2.0) * scale_, scale_, T(2.0) * scale_};
        scales[0] = hep::scales<T>{T(0.5) * scale_, scale_,          scale_};
        scales[0] = hep::scales<T>{         scale_, scale_, T(0.5) * scale_};
        scales[0] = hep::scales<T>{T(2.0) * scale_, scale_,          scale_};
        scales[0] = hep::scales<T>{         scale_, scale_, T(2.0) * scale_};
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

template <typename T>
std::size_t static_scale_function<T>::count() const
{
    switch (variation_)
    {
    case scale_variation::single_scale: return 1;
    case scale_variation::three_point:  return 3;
    case scale_variation::seven_point:  return 7;

    default:
        assert( false );
    }
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class static_scale_function<double>;

}
