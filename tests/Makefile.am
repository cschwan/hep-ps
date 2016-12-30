AM_CPPFLAGS = -I$(top_srcdir)/include -I$(srcdir) -DHEP_TYPE_T=double
AM_DEFAULT_SOURCE_EXT = .cpp

libcatch_la_SOURCES = catch.cpp
check_LTLIBRARIES = libcatch.la

noinst_HEADERS = catch.hpp test_structures.hpp

check_PROGRAMS = \
	test_fortran_helper \
	test_initial_state \
	test_lusifer_phase_space_generator \
	test_observables_born_like \
	test_observables_real

LDADD = $(check_LTLIBRARIES) $(top_builddir)/src/libhep-ps.la

TESTS = $(check_PROGRAMS)
