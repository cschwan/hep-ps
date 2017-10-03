AM_CPPFLAGS = -I$(top_srcdir)/include -I$(srcdir) -DHEP_TYPE_T=double
AM_DEFAULT_SOURCE_EXT = .cpp

libcatch_la_SOURCES = catch.cpp
check_LTLIBRARIES = libcatch.la

noinst_HEADERS = catch.hpp test_structures.hpp

check_PROGRAMS = \
	test_born_integrand \
	test_cofferaa_phase_space_generator \
	test_fini_integrand \
	test_fortran_helper \
	test_headers \
	test_lusifer_diagram_generator \
	test_lusifer_phase_space_generator \
	test_lusifer_phase_space_generator2 \
	test_lusifer_phase_space_generator3 \
	test_permutation \
	test_phase_space_quality \
	test_rambo_phase_space_generator \
	test_real_integrand \
	test_trivial_cutter_and_recombiner

if ADDITIONAL_TESTS
check_PROGRAMS += \
	test_p_type_jet_algorithm
endif

LDADD = $(check_LTLIBRARIES) $(top_builddir)/src/libhep-ps.la

TESTS = $(check_PROGRAMS)
