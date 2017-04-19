#include "hep/ps/lusifer_phase_space_generator.hpp"
#include "hep/ps/p_type_jet_algorithm.hpp"
#include "hep/ps/random_numbers.hpp"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

#include "catch.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

using T = double;

hep::lusifer_constants<T> constants(
	T(125.09), T(4.0e-3),
	T(174.2), T(1.41),
	T(80.385), T(2.085),
	T(91.1876), T(2.4952)
);

TEST_CASE("comparison against FastJet", "[p_type_jet_algorithm]")
{
	hep::lusifer_phase_space_generator<T> psg(
		"sq~uq ne el~nm mu~dq cq~gl ",
		constants
	);

	T const R = T(0.7);
	T const cmf_energy = T(1000.0);

	std::vector<std::size_t> proto_jets = { 6, 7, 8 };

	std::vector<T> numbers(psg.dimensions());
	std::vector<T> momenta(psg.map_dimensions());
	std::vector<T> recombined_momenta;
	std::mt19937 rng;

	std::vector<fastjet::PseudoJet> fj_momenta;
	std::vector<fastjet::PseudoJet> fj_recombined_momenta;
	std::vector<fastjet::PseudoJet> ps_recombined_momenta;

	for (double p : { -1.0, 0.0, 1.0 })
	{
		hep::p_type_jet_algorithm<T> jet_algorithm(T(p), R);

		fastjet::JetDefinition jet_definition(fastjet::genkt_algorithm, R, p);

		for (std::size_t channel = 0; channel != psg.channels(); ++channel)
		{
			for (std::size_t i = 0; i != 100; ++i)
			{
				std::generate(numbers.begin(), numbers.end(), [&]() {
					return std::generate_canonical<T,
						std::numeric_limits<T>::max_digits10>(rng);
				});

				hep::random_numbers<T> random_numbers(numbers);
				psg.generate(random_numbers, momenta, cmf_energy, channel);

				auto const recombinations = jet_algorithm.recombine(
					momenta,
					recombined_momenta,
					proto_jets,
					proto_jets.size()
				);

				fj_momenta.clear();

				for (auto const j : proto_jets)
				{
					fj_momenta.emplace_back(
						momenta.at(4 * j + 1),
						momenta.at(4 * j + 2),
						momenta.at(4 * j + 3),
						momenta.at(4 * j + 0)
					);
				}

				fj_recombined_momenta = fastjet::sorted_by_pt(
					fastjet::ClusterSequence(fj_momenta, jet_definition).
					inclusive_jets());

				// check if the number of recombined momenta agrees with FastJet
				CHECK( (proto_jets.size() - recombinations) ==
					fj_recombined_momenta.size() );

				ps_recombined_momenta.clear();

				for (std::size_t j = proto_jets.front();
					j != recombined_momenta.size() / 4; ++j)
				{
					ps_recombined_momenta.emplace_back(
						recombined_momenta.at(4 * j + 1),
						recombined_momenta.at(4 * j + 2),
						recombined_momenta.at(4 * j + 3),
						recombined_momenta.at(4 * j + 0)
					);
				}

				ps_recombined_momenta =
					fastjet::sorted_by_pt(ps_recombined_momenta);

				for (std::size_t j = 0; j != ps_recombined_momenta.size(); ++j)
				{
					CHECK( ps_recombined_momenta.at(j).e()  ==
						fj_recombined_momenta.at(j).e() );
					CHECK( ps_recombined_momenta.at(j).px() ==
						fj_recombined_momenta.at(j).px() );
					CHECK( ps_recombined_momenta.at(j).py() ==
						fj_recombined_momenta.at(j).py() );
					CHECK( ps_recombined_momenta.at(j).pz() ==
						fj_recombined_momenta.at(j).pz() );
				}
			}
		}
	}
}
