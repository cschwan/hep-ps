add_project_arguments('-DHAVE_CONFIG_H', language : 'cpp')

conf = configuration_data()
# get rid of this and use Fortran's C interfaces instead
conf.set('F77_FUNC_(name,NAME)', 'name ## _')

if get_option('openloops')
	conf.set('HAVE_OPENLOOPS', '1')
endif

configure_file(output : 'config.hpp', configuration : conf)

deps = []

if get_option('openloops')
	deps += [ meson.get_compiler('cpp').find_library('openloops') ]
endif

sources = [
	'boost.cpp',
	'cofferaa_generator.f',
	'cofferaa_extra.f',
	'cofferaa_interfaces.hpp',
	'cofferaa_phase_space_generator.cpp',
	'cs_subtraction.cpp',
	'diagram.cpp',
	'find_hep_ps.cpp',
	'fortran_helper.cpp',
	'generator.f',
	'hadron_hadron_psg_adapter.hpp',
	'initial_state.cpp',
	'inv_idx.cpp',
	'kaellen.cpp',
	'list_phase_space_generator.cpp',
	'lusifer_constants.cpp',
	'lusifer_diagram_generator.cpp',
	'lusifer_extra.f',
	'lusifer_interfaces.hpp',
	'lusifer_phase_space_generator.cpp',
	'ol_interface.cpp',
	'p_type_jet_algorithm.cpp',
	'p_type_photon_parton_recombiner.cpp',
	'parton.cpp',
	'parton_dfs.cpp',
	'parton_df_unc.cpp',
	'pdg_functions.cpp',
	'phase_space_generator.cpp',
	'phase_space_point.cpp',
	'propagator.cpp',
	'rambo_phase_space_generator.cpp',
	'static_scale_function.cpp',
	'silence_third_party_libs.cpp',
	'trivial_cutter.cpp',
	'trivial_recombiner.cpp',
	'vertex.cpp'
]

deps += dependency('gsl')
deps += dependency('lhapdf', version : '>=6')

# remove this line once the following bug is fixed:
# https://github.com/mesonbuild/meson/issues/2685
deps += meson.get_compiler('cpp').find_library('gfortran')

libhep_ps = library(
	'hep-ps',
	sources,
	dependencies : deps,
	include_directories : incdir,
	install : true
)
