#include "hep/ps/initial_state.hpp"
#include "hep/ps/pp_luminosities.hpp"

#include <LHAPDF/PDF.h>

#include <cassert>
#include <utility>

namespace hep
{

template <typename T>
class pp_luminosities<T>::impl
{
public:
	impl(std::string const& pdf_name, std::size_t pdf_member);

	std::unique_ptr<LHAPDF::PDF> pdf;
	std::vector<double> xfx1;
	std::vector<double> xfx2;
};

template <typename T>
pp_luminosities<T>::impl::impl(
	std::string const& pdf_name,
	std::size_t pdf_member
)
	: pdf(LHAPDF::mkPDF(pdf_name, pdf_member))
	, xfx1(13)
	, xfx2(13)
{
}

template <typename T>
pp_luminosities<T>::pp_luminosities(
	std::string const& pdf_name,
	std::size_t pdf_member
)
	// TODO: replace with `make_unique` in C++14
//	: pimpl(std::make_unique<impl>(pdf_name, pdf_member))
	: pimpl(new impl(pdf_name, pdf_member))
{
}

template <typename T>
pp_luminosities<T>::pp_luminosities(pp_luminosities&& luminosities)
	: pimpl(std::move(luminosities.pimpl))
{
}

template <typename T>
pp_luminosities<T>::~pp_luminosities() = default;

template <typename T>
T pp_luminosities<T>::alphas(T scale)
{
	return pimpl->pdf->alphasQ(scale);
}

template <typename T>
initial_state_array<T> pp_luminosities<T>::pdfs(T x1, T x2, T scale)
{
	initial_state_array<T> result;

	auto& xfx1 = pimpl->xfx1;
	auto& xfx2 = pimpl->xfx2;

	pimpl->pdf->xfxQ(static_cast <double> (x1), scale, xfx1);
	pimpl->pdf->xfxQ(static_cast <double> (x2), scale, xfx2);

	// TODO: is it possible to get these values from LHAPDF?
	constexpr std::size_t strange = -3 + 6; // anti strange
	constexpr std::size_t down    = -1 + 6; // anti down
	constexpr std::size_t gluon   =  0 + 6;
	constexpr std::size_t up      =  2 + 6;
	constexpr std::size_t charm   =  4 + 6;

	assert( strange < xfx1.size() );
	assert( down    < xfx1.size() );
	assert( gluon   < xfx1.size() );
	assert( up      < xfx1.size() );
	assert( charm   < xfx1.size() );

	for (auto& xfx : xfx1)
	{
		xfx /= x1;
	}

	for (auto& xfx : xfx2)
	{
		xfx /= x2;
	}

	result.set(initial_state::q43_uu, T(0.5 * xfx1[up]      * xfx2[up]));
	result.set(initial_state::q43_cc, T(0.5 * xfx1[charm]   * xfx2[charm]));
	result.set(initial_state::q23_dd, T(0.5 * xfx1[down]    * xfx2[down]));
	result.set(initial_state::q23_ss, T(0.5 * xfx1[strange] * xfx2[strange]));

	result.set(initial_state::x43_uu, result.get(initial_state::q43_uu));
	result.set(initial_state::x43_cc, result.get(initial_state::q43_cc));
	result.set(initial_state::x23_dd, result.get(initial_state::q23_dd));
	result.set(initial_state::x23_ss, result.get(initial_state::q23_ss));

	// Q=4/3
	//

	result.set(initial_state::q43_uc, T(xfx1[up]      * xfx2[charm]));
	result.set(initial_state::q43_cu, T(xfx1[charm]   * xfx2[up]));

	// Q=3/3
	//

	result.set(initial_state::q33_ud, T(xfx1[up]      * xfx2[down]));
	result.set(initial_state::q33_du, T(xfx1[down]    * xfx2[up]));
	result.set(initial_state::q33_us, T(xfx1[up]      * xfx2[strange]));
	result.set(initial_state::q33_su, T(xfx1[strange] * xfx2[up]));
	result.set(initial_state::q33_cd, T(xfx1[charm]   * xfx2[down]));
	result.set(initial_state::q33_dc, T(xfx1[down]    * xfx2[charm]));
	result.set(initial_state::q33_cs, T(xfx1[charm]   * xfx2[strange]));
	result.set(initial_state::q33_sc, T(xfx1[strange] * xfx2[charm]));

	// Q=2/3
	//

	result.set(initial_state::q23_ds, T(xfx1[down]    * xfx2[strange]));
	result.set(initial_state::q23_sd, T(xfx1[strange] * xfx2[down]));

	// Q=1/3
	//

	result.set(initial_state::q13_dg, T(xfx1[down]    * xfx2[gluon]));
	result.set(initial_state::q13_gd, T(xfx1[gluon]   * xfx2[down]));
	result.set(initial_state::q13_sg, T(xfx1[strange] * xfx2[gluon]));
	result.set(initial_state::q13_gs, T(xfx1[gluon]   * xfx2[strange]));

	// Q=2/3
	//

	result.set(initial_state::q23_ug, T(xfx1[up]      * xfx2[gluon]));
	result.set(initial_state::q23_gu, T(xfx1[gluon]   * xfx2[up]));
	result.set(initial_state::q23_cg, T(xfx1[charm]   * xfx2[gluon]));
	result.set(initial_state::q23_gc, T(xfx1[gluon]   * xfx2[charm]));

	// TODO: add remaining luminosities/think of how to generalize it

	return result;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class pp_luminosities<double>;

}
