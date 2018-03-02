#include "hep/ps/list_phase_space_generator.hpp"

#include <cassert>
#include <cstddef>

namespace
{

template <typename T>
class list_psg : public hep::phase_space_generator<T>
{
public:
	list_psg(std::vector<hep::list_phase_space_point<T>> const& list);

	std::size_t channels() const override;

	T densities(std::vector<T>& densities) override;

	std::size_t dimensions() const override;

	void generate(
		std::vector<T> const&,
		std::vector<T>& momenta,
		std::size_t channel
	) override;

	hep::luminosity_info<T> info() const override;

	std::size_t map_dimensions() const override;

private:
	std::size_t index_;
	std::vector<hep::list_phase_space_point<T>> list_;
};

template <typename T>
list_psg<T>::list_psg(std::vector<hep::list_phase_space_point<T>> const& list)
	: index_{0}
	, list_(list)
{
}

template <typename T>
std::size_t list_psg<T>::channels() const
{
	return 1;
}

template <typename T>
T list_psg<T>::densities(std::vector<T>& densities)
{
	densities.assign(1, T(1.0));

	return T(1.0);
}

template <typename T>
std::size_t list_psg<T>::dimensions() const
{
	return 0;
}

template <typename T>
void list_psg<T>::generate(
	std::vector<T> const&,
	std::vector<T>& momenta,
	std::size_t channel
) {
	assert( channel == 0 );

	auto const& phase_space = list_.at(index_).phase_space();

	momenta.assign(phase_space.begin(), phase_space.end());
	index_ = (index_ + 1) % list_.size();
}

template <typename T>
hep::luminosity_info<T> list_psg<T>::info() const
{
	if (index_ != 0)
	{
		return list_.at(index_ - 1).info();
	}
	else
	{
		return list_.back().info();
	}
}

template <typename T>
std::size_t list_psg<T>::map_dimensions() const
{
	return list_.front().phase_space().size();
}

template class list_psg<double>;

}

namespace hep
{

template <typename T>
list_phase_space_point<T>::list_phase_space_point(
	luminosity_info<T> const& info,
	std::vector<T> const& phase_space
)
	: info_{info}
	, phase_space_(phase_space)
{
}

template <typename T>
luminosity_info<T> const& list_phase_space_point<T>::info() const
{
	return info_;
}

template <typename T>
std::vector<T> const& list_phase_space_point<T>::phase_space() const
{
	return phase_space_;
}

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_list_phase_space_generator(
	std::vector<list_phase_space_point<T>> const& list_of_phase_space_points
) {
	return std::make_unique<list_psg<T>>(list_of_phase_space_points);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class list_phase_space_point<double>;

template std::unique_ptr<phase_space_generator<double>>
make_list_phase_space_generator(
	std::vector<list_phase_space_point<double>> const&
);

}
