project(
	'hep-ps',
	'cpp',
	'fortran',
	version : '0.0.1',
	license : 'GPL3',
	default_options : 'cpp_std=c++14',
	meson_version : '>=0.42'
)

hep_mc = dependency(
	'hep-mc',
	fallback : ['hep-mc', 'hep_mc_dep'],
	version : '>=0.7'
)

subdir('include')

if get_option('doxygen')
	subdir('doc')
endif

subdir('src')
subdir('tests')

pkg = import('pkgconfig')
pkg.generate(
	name : meson.project_name(),
	description : 'A C++ Library of Phase Space Integrands for High Energy Physics',
	libraries : libhep_ps,
	version : meson.project_version()
)

hep_ps_dep = declare_dependency(
	dependencies : hep_mc,
	include_directories : incdir,
	link_with : libhep_ps
)
