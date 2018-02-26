#include "hep/ps/list_phase_space_generator.hpp"

#include <cassert>
#include <cstddef>

namespace
{

template <typename T>
class list_psg : public hep::phase_space_generator<T>
{
public:
	list_psg(std::vector<std::vector<T>> const& list);

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
	std::vector<std::vector<T>> list_;
};

template <typename T>
list_psg<T>::list_psg(std::vector<std::vector<T>> const& list)
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

	momenta.assign(list_.at(index_).begin(), list_.at(index_).end());
	index_ = (index_ + 1) % list_.size();
}

template <typename T>
hep::luminosity_info<T> list_psg<T>::info() const
{
	return hep::luminosity_info<T>{};
}

template <typename T>
std::size_t list_psg<T>::map_dimensions() const
{
	return list_.front().size();
}

template class list_psg<double>;

}

namespace hep
{

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_list_phase_space_generator(
	std::vector<std::vector<T>> const& list_of_phase_space_points
) {
	return std::make_unique<list_psg<T>>(list_of_phase_space_points);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template std::unique_ptr<phase_space_generator<double>>
make_list_phase_space_generator(std::vector<std::vector<double>> const&);

}
