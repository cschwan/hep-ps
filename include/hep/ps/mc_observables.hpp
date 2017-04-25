#ifndef HEP_PS_MC_OBSERVABLES_HPP
#define HEP_PS_MC_OBSERVABLES_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2017  Christopher Schwan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hep/ps/initial_state_set.hpp"
#include "hep/ps/observables.hpp"
#include "hep/ps/observables_born.hpp"
#include "hep/ps/observables_fini.hpp"
#include "hep/ps/observables_real.hpp"
#include "hep/ps/mc_distributions.hpp"

#include "hep/mc/multi_channel_point.hpp"
#include "hep/mc/projector.hpp"

#include <memory>
#include <utility>

namespace
{

template <typename T, typename P>
hep::luminosity_info<T> info(hep::multi_channel_point2<T, P> const& x)
{
	return x.map().info();
}

template <typename T, typename P>
hep::luminosity_info<T> info(
	hep::multi_channel_point2<T, std::reference_wrapper<P>> const& x
) {
	return x.map().get().info();
}

}

namespace hep
{

template <typename T>
class mc_observables
{
public:
	mc_observables(
		std::unique_ptr<observables<T>>&& observables,
		initial_state_set set
	)
		: observables_(std::move(observables))
		, distributions_(nullptr)
		, set_(set)
	{
	}

	template <typename P>
	T operator()(
		hep::multi_channel_point2<T, P> const& point,
		hep::projector<T>& projector
	) {
		if (distributions_ == nullptr)
		{
			distributions_ = dynamic_cast <mc_distributions<T>*>
				(&observables_->distributions());
		}

		distributions_->set_projector(projector);
		return observables_->eval(point.coordinates(), info(point), set_);
	}

private:
	std::unique_ptr<observables<T>> observables_;
	mc_distributions<T>* distributions_;
	initial_state_set set_;
};

template <class T, class M, class C, class R, class P, class S, class D>
inline std::unique_ptr<mc_observables<T>> make_mc_observables_born(
	M&& matrix_elements,
	C&& cuts,
	R&& recombiner,
	P&& pdf,
	S&& scale_setter,
	D&& distributions,
	T hbarc2,
	initial_state_set set
) {
	return std::unique_ptr<mc_observables<T>>(new mc_observables<T>(
		make_observables_born(
			std::forward<M>(matrix_elements),
			std::forward<C>(cuts),
			std::forward<R>(recombiner),
			std::forward<P>(pdf),
			std::forward<S>(scale_setter),
			std::forward<D>(distributions),
			hbarc2
		),
		set
	));
}

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
inline std::unique_ptr<mc_observables<T>> make_mc_observables_fini(
	M&& matrix_elements,
	S&& subtraction,
	C&& cuts,
	R&& recombiner,
	P&& pdf,
	U&& scale_setter,
	D&& distributions,
	T hbarc2,
	bool insertion2,
	initial_state_set set
) {
	return std::unique_ptr<mc_observables<T>>(new mc_observables<T>(
		make_observables_fini(
			std::forward<M>(matrix_elements),
			std::forward<S>(subtraction),
			std::forward<C>(cuts),
			std::forward<R>(recombiner),
			std::forward<P>(pdf),
			std::forward<U>(scale_setter),
			std::forward<D>(distributions),
			hbarc2,
			insertion2
		),
		set
	));
}

template <class T, class M, class S, class C, class R, class P, class U,
	class D>
inline std::unique_ptr<mc_observables<T>> make_mc_observables_real(
	M&& matrix_elements,
	S&& subtraction,
	C&& cuts,
	R&& recombiner,
	P&& pdf,
	U&& scale_setter,
	D&& distributions,
	T hbarc2,
	T alpha_min,
	initial_state_set set
) {
	auto observable = make_observables_real(
		std::forward<M>(matrix_elements),
		std::forward<S>(subtraction),
		std::forward<C>(cuts),
		std::forward<R>(recombiner),
		std::forward<P>(pdf),
		std::forward<U>(scale_setter),
		std::forward<D>(distributions),
		hbarc2,
		alpha_min
	);

	return std::unique_ptr<mc_observables<T>>(new mc_observables<T>(
		std::move(observable),
		set
	));
}

}

#endif
