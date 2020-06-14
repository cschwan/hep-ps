project(
    'hep-ps',
    'cpp',
    'fortran',
    version : '0.0.1',
    license : 'GPL3',
    default_options : 'cpp_std=c++14',
    meson_version : '>=0.51'
)

hep_mc_dep = dependency(
    'hep-mc',
    fallback : ['hep-mc', 'hep_mc_dep'],
    version : '>=0.7'
)

span_lite_dep = dependency(
    'span-lite',
    method : 'cmake',
    modules : 'nonstd::span-lite',
    fallback : ['nonstd-span-lite', 'nonstd_span_lite_dep'],
    version : '>=0.5.0'
)

subdir('include')

if get_option('doxygen')
    subdir('doc')
endif

subdir('src')

pkg = import('pkgconfig')
pkg.generate(
    name : meson.project_name(),
    description : 'A C++ Library of Phase Space Integrands for High Energy Physics',
    libraries : libhep_ps,
    version : meson.project_version()
)

hep_ps_dep = declare_dependency(
    dependencies : hep_mc_dep,
    include_directories : incdir,
    link_with : libhep_ps
)

# tests depend on `hep_ps_dep`
subdir('tests')
