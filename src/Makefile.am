AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/include

lib_LTLIBRARIES = libhep-ps.la
libhep_ps_la_SOURCES = \
	boost.cpp \
	cofferaa_generator.f \
	cofferaa_extra.f \
	cofferaa_interfaces.hpp \
	cofferaa_phase_space_generator.cpp \
	cs_subtraction.cpp \
	diagram.cpp \
	ew_cs_subtraction.cpp \
	find_hep_ps.cpp \
	fortran_helper.cpp \
	generator.f \
	hadron_hadron_psg_adapter.hpp \
	inv_idx.cpp \
	lusifer_constants.cpp \
	lusifer_diagram_generator.cpp \
	lusifer_extra.f \
	lusifer_interfaces.hpp \
	lusifer_phase_space_generator.cpp \
	p_type_jet_algorithm.cpp \
	p_type_photon_parton_recombiner.cpp \
	parton.cpp \
	parton_dfs.cpp \
	parton_df_unc.cpp \
	phase_space_generator.cpp \
	phase_space_point.cpp \
	propagator.cpp \
	rambo_phase_space_generator.cpp \
	silence_third_party_libs.cpp \
	trivial_cutter.cpp \
	trivial_recombiner.cpp \
	vertex.cpp

libhep_ps_la_LIBADD = $(FLIBS)
