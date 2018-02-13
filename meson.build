project(
	'hep-ps',
	'cpp',
	'fortran',
	version : '0.0.1',
	license : 'GPL-3',
	default_options : 'cpp_std=c++14',
	meson_version : '>=0.42'
)

hep_mc = dependency('hep-mc', fallback : ['hep-mc', 'hep_mc_dep'])

if get_option('doxygen')
	subdir('doc')
endif

subdir('include')
subdir('src')
subdir('tests')

hep_ps_dep = declare_dependency(
	dependencies : hep_mc,
	include_directories : incdir,
	link_with : libhep_ps
)
