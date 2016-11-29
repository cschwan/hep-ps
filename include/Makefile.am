nobase_include_HEADERS = \
	hep/ps.hpp \
	hep/ps/constants.hpp \
	hep/ps/cut_result.hpp \
	hep/ps/dipole.hpp \
	hep/ps/dipole_type.hpp \
	hep/ps/p_type_jet_algorithm.hpp \
	hep/ps/particle_type.hpp \
	hep/ps/scales.hpp

# create a timestamp file so we can depend on it instead of depending on every
# single header
stamp-headers: $(nobase_include_HEADERS)
	echo "timestamp for headers" > $@

# build this by default
all-local: stamp-headers

MOSTLYCLEANFILES = stamp-headers
