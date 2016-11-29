AM_CPPFLAGS = -I$(top_srcdir)/include -I$(srcdir) -DHEP_TYPE_T=double
AM_DEFAULT_SOURCE_EXT = .cpp

libcatch_la_SOURCES = catch.cpp
check_LTLIBRARIES = libcatch.la

noinst_HEADERS = catch.hpp test_structures.hpp

check_PROGRAMS = \
	test_real_minus_dipoles

LDADD = $(check_LTLIBRARIES) $(top_builddir)/src/libhep-ps.la

TESTS = $(check_PROGRAMS)
