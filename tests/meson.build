libcatch = static_library(
	'catch',
	'catch.hpp',
	'catch.cpp'
)

test_born_integrand = executable(
	'test_born_integrand',
	'test_born_integrand.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : [ hep_mc ],
	link_with : [ libcatch, libhep_ps ]
)

test_cofferaa_phase_space_generator = executable(
	'test_cofferaa_phase_space_generator',
	'test_cofferaa_phase_space_generator.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_fini_integrand = executable(
	'test_fini_integrand',
	'test_fini_integrand.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_fortran_helper = executable(
	'test_fortran_helper',
	'test_fortran_helper.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	link_with : [ libcatch, libhep_ps ]
)

test_headers = executable(
	'test_headers',
	'test_headers.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_list_phase_space_generator = executable(
	'test_list_phase_space_generator',
	'test_list_phase_space_generator.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_lusifer_diagram_generator = executable(
	'test_lusifer_diagram_generator',
	'test_lusifer_diagram_generator.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	link_with : [ libcatch, libhep_ps ]
)

test_lusifer_phase_space_generator = executable(
	'test_lusifer_phase_space_generator',
	'test_lusifer_phase_space_generator.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_lusifer_phase_space_generator2 = executable(
	'test_lusifer_phase_space_generator2',
	'test_lusifer_phase_space_generator2.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_lusifer_phase_space_generator3 = executable(
	'test_lusifer_phase_space_generator3',
	'test_lusifer_phase_space_generator3.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_parton = executable(
	'test_parton',
	'test_parton.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	link_with : [ libcatch, libhep_ps ]
)

test_permutation = executable(
	'test_permutation',
	'test_permutation.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	link_with : [ libcatch, libhep_ps ]
)

test_phase_space_quality = executable(
	'test_phase_space_quality',
	'test_phase_space_quality.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_rambo_phase_space_generator = executable(
	'test_rambo_phase_space_generator',
	'test_rambo_phase_space_generator.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_real_integrand = executable(
	'test_real_integrand',
	'test_real_integrand.cpp',
	include_directories : incdir,
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	link_with : [ libcatch, libhep_ps ]
)

test_trivial_cutter_and_recombiner = executable(
	'test_trivial_cutter_and_recombiner',
	'test_trivial_cutter_and_recombiner.cpp',
	cpp_args : '-DHEP_TYPE_T=double',
	dependencies : hep_mc,
	include_directories : incdir,
	link_with : [ libcatch, libhep_ps ]
)

test('test_born_integrand', test_born_integrand)
test('test_cofferaa_phase_space_generator', test_cofferaa_phase_space_generator)
test('test_fini_integrand', test_fini_integrand)
test('test_fortran_helper', test_fortran_helper)
test('test_headers', test_headers)
test('test_list_phase_space_generator', test_list_phase_space_generator)
test('test_lusifer_diagram_generator', test_lusifer_diagram_generator)
test('test_lusifer_phase_space_generator', test_lusifer_phase_space_generator)
test('test_lusifer_phase_space_generator2', test_lusifer_phase_space_generator2)
test('test_lusifer_phase_space_generator3', test_lusifer_phase_space_generator3)
test('test_parton', test_parton)
test('test_permutation', test_permutation)
test('test_phase_space_quality', test_phase_space_quality, timeout : 120)
test('test_rambo_phase_space_generator', test_rambo_phase_space_generator)
test('test_real_integrand', test_real_integrand)
test('test_trivial_cutter_and_recombiner', test_trivial_cutter_and_recombiner)

if get_option('additional-tests')
	fastjet = dependency('fastjet')

	test_p_type_jet_algorithm = executable(
		'test_p_type_jet_algorithm',
		'test_p_type_jet_algorithm.cpp',
		cpp_args : '-DHEP_TYPE_T=double',
		dependencies : fastjet,
		include_directories : incdir,
		link_with : [ libcatch, libhep_ps ]
	)

	test('test_p_type_jet_algorithm', test_p_type_jet_algorithm)
endif
