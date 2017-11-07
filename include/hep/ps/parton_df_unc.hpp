#ifndef HEP_PS_PARTON_DF_UNC_HPP
#define HEP_PS_PARTON_DF_UNC_HPP

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace hep
{

/// Structure that captures the symmetric and asymmetric PDF uncertainties.
template <typename T>
struct pdf_unc
{
	/// Symmetric PDF uncertainties.
	T sym;

	/// Asymmetric PDF uncertainties (lower value).
	T neg;

	/// Asymmetric PDF uncertainties (upper value).
	T pos;
};

/// Class that calculates the PDF uncertainty for a given PDF set.
template <typename T>
class parton_df_unc
{
public:
	/// Constructor.
	parton_df_unc(std::string const& name, T cl = T());

	/// Move constructor.
	parton_df_unc(parton_df_unc&& pdf_unc);

	/// Destructor.
	~parton_df_unc();

	/// Calculates the PDF uncertainty for the given `values`. The first element
	/// of `values` must be the central prediction, all other values the ones
	/// from the Error PDF set.
	pdf_unc<T> uncertainty(std::vector<T> const& values) const;

private:
	class impl;
	std::unique_ptr<impl> pimpl;
};

}

#endif
