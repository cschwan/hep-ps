# AC_CONFIG_MACRO_DIR needs >=2.58
AC_PREREQ([2.58])

# initialize configure
AC_INIT([hep-ps], [0.0.1], [], [hep-ps], [])

# put autotools files into a separate directory
AC_CONFIG_AUX_DIR([build-aux])

# initialize automake
AM_INIT_AUTOMAKE([foreign -Wall dist-bzip2 no-dist-gzip nostdinc])

# put aclocal macros into m4
AC_CONFIG_MACRO_DIR([m4])

# C++ is our language
AC_LANG([C++])

# needed for compatibility with non-POSIX systems
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# check for C++11 conforming compiler (autoconf-archive macro)
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])

# we also need a Fortran 77 compiler
AC_PROG_F77

# needed for interoperability with C++
AC_F77_LIBRARY_LDFLAGS

# needed to find out FORTAN signatures
AC_F77_WRAPPERS

# initialize libtool
LT_INIT

# point to a source file in order to let ./configure find itself
AC_CONFIG_SRCDIR([src/find_hep_ps.cpp])

# don't pollute the command line with compiler defs, write them into a header
# that we dont need
AC_CONFIG_HEADERS([config.hpp])

# check for hep-mc (needed for some tests)
AC_CHECK_HEADERS([hep/mc.hpp])
AM_CONDITIONAL([HAVE_HEP_MC_HPP], [test -n "$ac_cv_header_hep_mc_hpp"])

# check for FastJet (needed for a test)
AC_CHECK_HEADERS([fastjet/ClusterSequence.hh],
	AC_SEARCH_LIBS([_ZN7fastjet15ClusterSequenceD1Ev], [fastjet]))
AM_CONDITIONAL([HAVE_FASTJET_CLUSTER_SEQUENCE_HH],
	[test -n "$ac_cv_header_fastjet_ClusterSequence_hh"])

# find LHAPDF header and library
AC_CHECK_HEADERS([LHAPDF/LHAPDF.h])
AC_SEARCH_LIBS([_ZN6LHAPDF5mkPDFERKSsi], [LHAPDF], [],
	AC_MSG_ERROR([LHAPDF is a required dependency]))

# if the documentation is enabled find doxygen binary
AC_ARG_ENABLE([doxygen], AS_HELP_STRING([--enable-doxygen],
	[Enable generation of Doxygen documentation]))
AC_ARG_VAR([DOXYGEN], [Doxygen documentation generator])
AS_IF([test "x$enable_doxygen" = "xyes"], [
	AC_PATH_PROGS([DOXYGEN], [doxygen])
	AS_IF([test -z "$DOXYGEN"], AC_MSG_ERROR([Doxygen not found but requested]))
])
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])

# configure Makefiles
AC_CONFIG_FILES([
	Makefile
	doc/Makefile
	include/Makefile
	src/Makefile
	tests/Makefile
])

# configure remaining files
AC_CONFIG_FILES([doc/Doxyfile])

# write configure script
AC_OUTPUT
