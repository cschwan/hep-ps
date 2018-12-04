catch_subproj = subproject('catch')
catch_dep = catch_subproj.get_variable('catch_dep')

# make it a static library because `catch.cpp` contains `main`
libcatch = static_library('catch', 'catch.cpp', dependencies : catch_dep)

libcatch_dep = declare_dependency(dependencies : catch_dep, link_with : libcatch)

programs = [
    'regression_test_lusifer_psg_ww_scattering',
    'regression_test_lusifer_psg_wz_scattering',
    'regression_test_lusifer_psg_www_production',
    'regression_test_lusifer_psg_wwzph_production',
    'regression_test_lusifer_psg_wzj_production',
    'test_ps_invariant_commutation',
    'test_born_integrand',
    'test_cofferaa_phase_space_generator',
    'test_dipole',
    'test_fini_integrand',
    'test_fortran_helper',
    'test_generate_dipole',
    'test_headers',
    'test_initial_state',
    'test_insertion_term',
    'test_list_phase_space_generator',
    'test_lusifer_phase_space_generator',
    'test_lusifer_ps_channels',
    'test_ol_born_matrix_elements',
    'test_pdg_functions',
    'test_permutation',
    'test_phase_space_point',
    'test_phase_space_quality',
    'test_rambo_phase_space_generator',
    'test_real_integrand',
    'test_trivial_cutter_and_recombiner'
]

foreach p : programs
    exe = executable(
        p,
        p + '.cpp',
        cpp_args : '-DHEP_TYPE_T=double',
        dependencies : [ libcatch_dep, hep_mc_dep, hep_ps_dep ]
    )

    test(p, exe, timeout : 60)
endforeach

if get_option('additional-tests')
    # TODO: this doesn't work, fix it
    fastjet_dep = dependency('fastjet')

    test_p_type_jet_algorithm = executable(
        'test_p_type_jet_algorithm',
        'test_p_type_jet_algorithm.cpp',
        cpp_args : '-DHEP_TYPE_T=double',
        dependencies : [ fastjet_dep, libcatch_dep, hep_mc_dep, hep_ps_dep ]
    )

    test('test_p_type_jet_algorithm', test_p_type_jet_algorithm)
endif
