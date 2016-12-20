ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = NEWS.rst README.rst

# if doxygen is enabled include doc subdir
if HAVE_DOXYGEN
    DOCDIR = doc
endif

# `doc` must be built after `include`
SUBDIRS = include $(DOCDIR) src tests
