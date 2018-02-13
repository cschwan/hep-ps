conf = configuration_data()
# get rid of this and use Fortran's C interfaces instead
conf.set('F77_FUNC_(name,NAME)', 'name ## _')
configure_file(output : 'config.hpp', configuration : conf)
add_project_arguments('-DHAVE_CONFIG_H', language : 'cpp')

sources = [
	'boost.cpp',
	'cofferaa_generator.f',
	'cofferaa_extra.f',
	'cofferaa_interfaces.hpp',
	'cofferaa_phase_space_generator.cpp',
	'cs_subtraction.cpp',
	'diagram.cpp',
	'ew_cs_subtraction.cpp',
	'find_hep_ps.cpp',
	'fortran_helper.cpp',
	'generator.f',
	'hadron_hadron_psg_adapter.hpp',
	'initial_state.cpp',
	'inv_idx.cpp',
	'lusifer_constants.cpp',
	'lusifer_diagram_generator.cpp',
	'lusifer_extra.f',
	'lusifer_interfaces.hpp',
	'lusifer_phase_space_generator.cpp',
	'p_type_jet_algorithm.cpp',
	'p_type_photon_parton_recombiner.cpp',
	'parton.cpp',
	'parton_dfs.cpp',
	'parton_df_unc.cpp',
	'phase_space_generator.cpp',
	'phase_space_point.cpp',
	'propagator.cpp',
	'rambo_phase_space_generator.cpp',
	'silence_third_party_libs.cpp',
	'trivial_cutter.cpp',
	'trivial_recombiner.cpp',
	'vertex.cpp'
]

gsl = dependency('gsl')
lhapdf = dependency('lhapdf', version : '>=6')

# remove this line once the following bug is fixed:
# https://github.com/mesonbuild/meson/issues/2685
fortran = meson.get_compiler('cpp').find_library('gfortran')

libhep_ps = library(
	'hep-ps',
	sources,
	dependencies : [ hep_mc, fortran, gsl, lhapdf ],
	include_directories : incdir
)
