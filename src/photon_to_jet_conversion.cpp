#include "hep/ps/photon_to_jet_conversion.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace hep
{

template <typename T>
photon_to_jet_conversion<T>::photon_to_jet_conversion()
{
}

template <typename T>
photon_to_jet_conversion<T>::photon_to_jet_conversion(
    T alpha,
    T delta_alpha_hadr,
    T conversion_scale,
    T nc,
    std::size_t nf
)
    : alpha_{alpha}
    , delta_alpha_hadr_{delta_alpha_hadr}
    , conversion_scale_{conversion_scale}
    , nc_{nc}
    , nf_{T(nf)}
{
    if (nf_ == 0)
    {
        throw std::invalid_argument("number of flavors must be larger than zero");
    }
}

template <typename T>
bool photon_to_jet_conversion<T>::active() const
{
    return nf_ != 0;
}

template <typename T>
T photon_to_jet_conversion<T>::eval(T regularization_scale, T non_perturbative_factor) const
{
    using std::acos;

    T const pi = acos(T(-1.0));

    T result = T(5.0) / T(3.0) + T(2.0) * log(regularization_scale / conversion_scale_);
    result /= T(-3.0) * pi;
    result *= nc_;
    result += non_perturbative_factor / alpha_ / nf_ * delta_alpha_hadr_;

    return result;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class photon_to_jet_conversion<double>;

}
