add_project_arguments('-DHAVE_CONFIG_H', language : 'cpp')

conf = configuration_data()
# get rid of this and use Fortran's C interfaces instead
conf.set('F77_FUNC_(name,NAME)', 'name ## _')

if get_option('openloops')
    conf.set('HAVE_OPENLOOPS', '1')
endif

if get_option('recola')
    conf.set('HAVE_RECOLA', '1')
endif

configure_file(output : 'config.hpp', configuration : conf)

deps = [ hep_mc_dep ]

sources = [
    'cofferaa_generator.f90',
    'cofferaa_extra.f90',
    'cofferaa_interfaces.hpp',
    'cofferaa_phase_space_generator.cpp',
    'cs_subtraction.cpp',
    'd_subtraction.cpp',
    'dipole.cpp',
    'dipole_vertex.cpp',
    'dipole_veto.cpp',
    'find_hep_ps.cpp',
    'fortran_helper.cpp',
    'generate_dipole.cpp',
    'generator.f90',
    'hadron_hadron_psg_adapter.hpp',
    'initial_state.cpp',
    'insertion_term.cpp',
    'int_dipole.cpp',
    'list_phase_space_generator.cpp',
    'lusifer_constants.cpp',
    'lusifer_extra.f90',
    'lusifer_interfaces.hpp',
    'lusifer_phase_space_generator.cpp',
    'lusifer_ps_channels.cpp',
    'lusifer_ps_functions.cpp',
    'nfpc.cpp',
    'ol_born_matrix_elements.cpp',
    'ol_int_dipoles.cpp',
    'ol_integrated_mes.cpp',
    'ol_interface.cpp',
    'ol_ioperator.cpp',
    'ol_real_matrix_elements.cpp',
    'p_type_jet_algorithm.cpp',
    'p_type_photon_parton_recombiner.cpp',
    'parton_dfs.cpp',
    'parton_df_unc.cpp',
    'pdg_functions.cpp',
    'phase_space_generator.cpp',
    'phase_space_point.cpp',
    'phase_space_tools.cpp',
    'photon_cone_recombiner.cpp',
    'photon_to_jet_conversion.cpp',
    'ps_functions.cpp',
    'psp.cpp',
    'rambo_phase_space_generator.cpp',
    'recombined_state.cpp',
    'static_scale_function.cpp',
    'suppress_banners.cpp',
    'trivial_cutter.cpp',
    'trivial_recombiner.cpp'
]

# TODO: make openloops and recola incompatible with each other
if get_option('recola')
    add_project_arguments(['-DGNU_NAME_MANGLING'], language: 'cpp')
    deps += [ meson.get_compiler('cpp').find_library('recola') ]
    sources += [ 'recola_openloops_bridge.cpp' ]
elif get_option('openloops')
    deps += [ meson.get_compiler('cpp').find_library('openloops') ]
endif

deps += dependency('gsl')
deps += dependency('lhapdf', version : '>=6')
deps += span_lite_dep

libhep_ps = library(
    'hep-ps',
    sources,
    dependencies : deps,
    include_directories : incdir,
    install : true
)
