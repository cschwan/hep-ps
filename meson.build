project(
	'hep-ps',
	'cpp',
	'fortran',
	version : '0.0.1',
	license : 'GPL-3',
	default_options : 'cpp_std=c++14',
	meson_version : '>=0.42'
)

add_project_arguments('-DHAVE_CONFIG_H', language : 'cpp')

#mpidep = dependency('mpi', language : 'cpp', required : true)
incdir = include_directories('include')

hep_mc = dependency('hep-mc', fallback : ['hep-mc', 'hep_mc_dep'])

if get_option('doxygen')
	subdir('doc')
endif

subdir('include')
subdir('src')
subdir('tests')
