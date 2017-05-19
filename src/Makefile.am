AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/include

lib_LTLIBRARIES = libhep-ps.la
libhep_ps_la_SOURCES = \
	boost.cpp \
	cofferaa_generator.f \
	cofferaa_extra.f \
	cofferaa_interfaces.hpp \
	cofferaa_phase_space_generator.cpp \
	cs_subtraction.cpp \
	find_hep_ps.cpp \
	fortran_helper.cpp \
	generator.f \
	hadron_hadron_psg_adapter.hpp \
	lusifer_constants.cpp \
	lusifer_extra.f \
	lusifer_interfaces.hpp \
	lusifer_phase_space_generator.cpp \
	p_type_jet_algorithm.cpp \
	proton_pdf.cpp \
	silence_third_party_libs.cpp \
	trivial_cutter.cpp \
	trivial_recombiner.cpp

libhep_ps_la_LIBADD = $(FLIBS)
