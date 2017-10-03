nobase_include_HEADERS = \
	hep/ps.hpp \
	hep/ps/abc_terms.hpp \
	hep/ps/alphas_calc.hpp \
	hep/ps/boost.hpp \
	hep/ps/cofferaa_phase_space_generator.hpp \
	hep/ps/constants.hpp \
	hep/ps/convolute.hpp \
	hep/ps/cs_subtraction.hpp \
	hep/ps/cut_result.hpp \
	hep/ps/diagram.hpp \
	hep/ps/dipole.hpp \
	hep/ps/dipole_invariants.hpp \
	hep/ps/dipole_type.hpp \
	hep/ps/dipole_with_set.hpp \
	hep/ps/enum.hpp \
	hep/ps/event_type.hpp \
	hep/ps/factorization_scheme.hpp \
	hep/ps/find_hep_ps.hpp \
	hep/ps/fortran_helper.hpp \
	hep/ps/initial_state.hpp \
	hep/ps/insertion_term.hpp \
	hep/ps/insertion_term_type.hpp \
	hep/ps/inv_idx.hpp \
	hep/ps/kaellen.hpp \
	hep/ps/luminosity_info.hpp \
	hep/ps/lusifer_constants.hpp \
	hep/ps/lusifer_diagram_generator.hpp \
	hep/ps/lusifer_phase_space_generator.hpp \
	hep/ps/neg_pos_results.hpp \
	hep/ps/non_zero_dipole.hpp \
	hep/ps/observables.hpp \
	hep/ps/observables_born.hpp \
	hep/ps/observables_fini.hpp \
	hep/ps/observables_real.hpp \
	hep/ps/p_type_jet_algorithm.hpp \
	hep/ps/particle_type.hpp \
	hep/ps/parton.hpp \
	hep/ps/parton_dfs.hpp \
	hep/ps/permutation.hpp \
	hep/ps/phase_space_generator.hpp \
	hep/ps/propagator.hpp \
	hep/ps/rambo_phase_space_generator.hpp \
	hep/ps/renormalization_scheme.hpp \
	hep/ps/scales.hpp \
	hep/ps/silence_third_party_libs.hpp \
	hep/ps/trivial_cutter.hpp \
	hep/ps/trivial_distributions.hpp \
	hep/ps/trivial_recombiner.hpp \
	hep/ps/vertex.hpp

# create a timestamp file so we can depend on it instead of depending on every
# single header
stamp-headers: $(nobase_include_HEADERS)
	echo "timestamp for headers" > $@

# build this by default
all-local: stamp-headers

MOSTLYCLEANFILES = stamp-headers
