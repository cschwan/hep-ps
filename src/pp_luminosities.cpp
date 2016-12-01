#include "hep/ps/initial_state.hpp"
#include "hep/ps/pp_luminosities.hpp"

#include <LHAPDF/PDF.h>

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

	pimpl->pdf->xfxQ(static_cast <double> (x1), scale, pimpl->xfx1);
	pimpl->pdf->xfxQ(static_cast <double> (x2), scale, pimpl->xfx2);

	// TODO: is it possible to get these values from LHAPDF?
	constexpr std::size_t strange      = -3 + 6; // anti strange
	constexpr std::size_t down         = -1 + 6; // anti down
	constexpr std::size_t gluon        =  0 + 6;
	constexpr std::size_t up           =  2 + 6;
	constexpr std::size_t charm        =  4 + 6;

	auto const& xfx1 = pimpl->xfx1;
	auto const& xfx2 = pimpl->xfx2;

	T const factor = T(1.0) / (x1 * x2);

	// Q=4/3
	//

	result.set(initial_state::q43_uc, factor * T(xfx1.at(up)      * xfx2.at(charm)));
	result.set(initial_state::q43_cu, factor * T(xfx1.at(charm)   * xfx2.at(up)));
	result.set(initial_state::q43_uu, factor * T(xfx1.at(up)      * xfx2.at(up)));
	result.set(initial_state::q43_cc, factor * T(xfx1.at(charm)   * xfx2.at(charm)));

	// Q=3/3
	//

	result.set(initial_state::q33_ud, factor * T(xfx1.at(up)      * xfx2.at(down)));
	result.set(initial_state::q33_du, factor * T(xfx1.at(down)    * xfx2.at(up)));
	result.set(initial_state::q33_us, factor * T(xfx1.at(up)      * xfx2.at(strange)));
	result.set(initial_state::q33_su, factor * T(xfx1.at(strange) * xfx2.at(up)));
	result.set(initial_state::q33_cd, factor * T(xfx1.at(charm)   * xfx2.at(down)));
	result.set(initial_state::q33_dc, factor * T(xfx1.at(down)    * xfx2.at(charm)));
	result.set(initial_state::q33_cs, factor * T(xfx1.at(charm)   * xfx2.at(strange)));
	result.set(initial_state::q33_sc, factor * T(xfx1.at(strange) * xfx2.at(charm)));

	// Q=2/3
	//

	result.set(initial_state::q23_ds, factor * T(xfx1.at(down)    * xfx2.at(strange)));
	result.set(initial_state::q23_sd, factor * T(xfx1.at(strange) * xfx2.at(down)));
	result.set(initial_state::q23_dd, factor * T(xfx1.at(down)    * xfx2.at(down)));
	result.set(initial_state::q23_ss, factor * T(xfx1.at(strange) * xfx2.at(strange)));

	// Q=1/3
	//

	result.set(initial_state::q13_dg, factor * T(xfx1.at(down)    * xfx2.at(gluon)));
	result.set(initial_state::q13_gd, factor * T(xfx1.at(gluon)   * xfx2.at(down)));
	result.set(initial_state::q13_sg, factor * T(xfx1.at(strange) * xfx2.at(gluon)));
	result.set(initial_state::q13_gs, factor * T(xfx1.at(gluon)   * xfx2.at(strange)));

	// Q=2/3
	//

	result.set(initial_state::q23_ug, factor * T(xfx1.at(up)      * xfx2.at(gluon)));
	result.set(initial_state::q23_gu, factor * T(xfx1.at(gluon)   * xfx2.at(up)));
	result.set(initial_state::q23_cg, factor * T(xfx1.at(charm)   * xfx2.at(gluon)));
	result.set(initial_state::q23_gc, factor * T(xfx1.at(gluon)   * xfx2.at(charm)));

	// TODO: add remaining luminosities/think of how to generalize it

	return result;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class pp_luminosities<double>;

}
