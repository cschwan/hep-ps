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

# initialize libtool
LT_INIT

# point to a source file in order to let ./configure find itself
AC_CONFIG_SRCDIR([src/find_hep_ps.cpp])

# don't pollute the command line with compiler defs, write them into a header
# that we dont need
AC_CONFIG_HEADERS([config.hpp])

# configure Makefiles and Doxyfile
AC_CONFIG_FILES([
	Makefile
	include/Makefile
	src/Makefile
])

# write configure script
AC_OUTPUT
