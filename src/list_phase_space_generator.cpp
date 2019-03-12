#include "hep/ps/list_phase_space_generator.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <istream>
#include <sstream>
#include <string>

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

    void generate(std::vector<T> const&, std::vector<T>& momenta, std::size_t channel) override;

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
    densities.assign(1, (index_ != 0) ? list_.at(index_ - 1).weight() : list_.back().weight());

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
    return (index_ != 0) ? list_.at(index_ - 1).info() : list_.back().info();
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
    std::vector<T> const& phase_space,
    T weight
)
    : info_{info}
    , phase_space_(phase_space)
    , weight_{weight}
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
T list_phase_space_point<T>::weight() const
{
    return weight_;
}

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_list_phase_space_generator(
    std::vector<list_phase_space_point<T>> const& list_of_phase_space_points
) {
    return std::make_unique<list_psg<T>>(list_of_phase_space_points);
}

template <typename T>
std::unique_ptr<phase_space_generator<T>> make_list_phase_space_generator(std::istream& stream)
{
    using std::log;

    std::string line;
    std::vector<T> numbers;
    std::vector<T> phase_space_point;
    std::size_t n = 0;
    std::size_t line_number = 0;
    std::vector<list_phase_space_point<T>> list;

    numbers.reserve(4);
    phase_space_point.reserve(40);

    bool read_point = true;

    while (std::getline(stream, line))
    {
        ++line_number;

        // skip empty and comment lines
        if ((line.size() == 0) || (line.at(0) == '#'))
        {
            continue;
        }

        std::stringstream tokens{line};

        do
        {
            T number;
            tokens >> number;
            numbers.push_back(number);
        }
        while (tokens.good());

        // if theres more than one number on a line, count it towards the phase space point
        if ((numbers.size() >= 4) && read_point)
        {
            phase_space_point.insert(phase_space_point.end(), numbers.begin(), numbers.end());
            numbers.clear();
        }
        else if (numbers.size() == 1)
        {
            if (!phase_space_point.empty())
            {
                if (n == 0)
                {
                    n = phase_space_point.size();
                }
                else if (n != phase_space_point.size())
                {
                    // TODO: format a proper error message and throw an exception
                    throw 0;
                }
            }

            read_point = false;
        }
        else if ((numbers.size() == 4) && !read_point)
        {
            T const x1 = numbers.at(0);
            T const x2 = numbers.at(1);
            T const energy = numbers.at(2);
            T const weight = numbers.at(3);
            T const rapidity_shift = T(0.5) * log(x1 / x2);

            list.emplace_back(luminosity_info<T>{x1, x2, energy * energy, rapidity_shift},
                phase_space_point, weight);
            phase_space_point.clear();
            numbers.clear();

            read_point = true;
        }
    }

    return make_list_phase_space_generator(list);
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template class list_phase_space_point<double>;

template std::unique_ptr<phase_space_generator<double>>
make_list_phase_space_generator(std::vector<list_phase_space_point<double>> const&);

template std::unique_ptr<phase_space_generator<double>>
make_list_phase_space_generator(std::istream& stream);

}
