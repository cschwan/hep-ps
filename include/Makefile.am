nobase_include_HEADERS = \
	hep/ps.hpp \
	hep/ps-mc.hpp \
	hep/ps/abc_terms.hpp \
	hep/ps/constants.hpp \
	hep/ps/cs_subtraction.hpp \
	hep/ps/cut_result.hpp \
	hep/ps/dipole.hpp \
	hep/ps/dipole_invariants.hpp \
	hep/ps/dipole_type.hpp \
	hep/ps/enum.hpp \
	hep/ps/event_type.hpp \
	hep/ps/factorization_scheme.hpp \
	hep/ps/fold.hpp \
	hep/ps/fortran_helper.hpp \
	hep/ps/hh_phase_space_generator.hpp \
	hep/ps/initial_state.hpp \
	hep/ps/initial_state_set.hpp \
	hep/ps/kaellen.hpp \
	hep/ps/luminosity_info.hpp \
	hep/ps/lusifer_phase_space_generator.hpp \
	hep/ps/mc_phase_space_adapter.hpp \
	hep/ps/neg_pos_results.hpp \
	hep/ps/observables_born_like.hpp \
	hep/ps/observables_real.hpp \
	hep/ps/p_type_jet_algorithm.hpp \
	hep/ps/permutation.hpp \
	hep/ps/pp_luminosities.hpp \
	hep/ps/particle_type.hpp \
	hep/ps/requires_cut.hpp \
	hep/ps/scales.hpp \
	hep/ps/silence_third_party_libs.hpp \
	hep/ps/trivial_cutter.hpp \
	hep/ps/trivial_distributions.hpp \
	hep/ps/trivial_recombiner.hpp

# create a timestamp file so we can depend on it instead of depending on every
# single header
stamp-headers: $(nobase_include_HEADERS)
	echo "timestamp for headers" > $@

# build this by default
all-local: stamp-headers

MOSTLYCLEANFILES = stamp-headers
