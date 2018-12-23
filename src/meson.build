add_project_arguments('-DHAVE_CONFIG_H', language : 'cpp')

conf = configuration_data()
# get rid of this and use Fortran's C interfaces instead
conf.set('F77_FUNC_(name,NAME)', 'name ## _')

if get_option('openloops')
    conf.set('HAVE_OPENLOOPS', '1')
endif

configure_file(output : 'config.hpp', configuration : conf)

deps = [ hep_mc_dep ]

if get_option('openloops')
    deps += [ meson.get_compiler('cpp').find_library('openloops') ]
endif

sources = [
    'cofferaa_generator.f90',
    'cofferaa_extra.f90',
    'cofferaa_interfaces.hpp',
    'cofferaa_phase_space_generator.cpp',
    'cs_subtraction.cpp',
    'dipole.cpp',
    'dipole_veto.cpp',
    'find_hep_ps.cpp',
    'fortran_helper.cpp',
    'generate_dipole.cpp',
    'generator.f90',
    'hadron_hadron_psg_adapter.hpp',
    'initial_state.cpp',
    'insertion_term.cpp',
    'list_phase_space_generator.cpp',
    'lusifer_constants.cpp',
    'lusifer_extra.f90',
    'lusifer_interfaces.hpp',
    'lusifer_phase_space_generator.cpp',
    'lusifer_ps_channels.cpp',
    'lusifer_ps_functions.cpp',
    'ol_born_matrix_elements.cpp',
    'ol_integrated_matrix_elements.cpp',
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
    'ps_functions.cpp',
    'psp.cpp',
    'rambo_phase_space_generator.cpp',
    'recombined_state.cpp',
    'static_scale_function.cpp',
    'suppress_banners.cpp',
    'trivial_cutter.cpp',
    'trivial_recombiner.cpp'
]

deps += dependency('gsl')
deps += dependency('lhapdf', version : '>=6')

libhep_ps = library(
    'hep-ps',
    sources,
    dependencies : deps,
    include_directories : incdir,
    install : true
)
