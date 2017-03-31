#ifdef HAVE_CONFIG_H
#include "config.hpp"
#endif

#define lusifer_initphasespace F77_FUNC(lusifer_initphasespace, LUSIFER_INITPHASESPACE)
#define lusifer_phasespace F77_FUNC(lusifer_phasespace, LUSIFER_PHASESPACE)
#define lusifer_density F77_FUNC(lusifer_density, LUSIFER_DENSITY)

#define lusifer_extra_max F77_FUNC(lusifer_extra_max, LUSIFER_EXTRA_MAX)
#define lusifer_extra_data F77_FUNC(lusifer_extra_data, LUSIFER_EXTRA_DATA)
#define lusifer_extra_set F77_FUNC(lusifer_extra_set, LUSIFER_EXTRA_SET)

extern "C"
{

void lusifer_initphasespace(
    char const* name,
    int* generator,
    int* lightfermions,
    int* includecuts,
    int* sout,
    int name_length
);

void lusifer_phasespace(
    double const* random,
    double* kbeam,
    double* k,
    double* x1,
    double* x2,
    double* g,
    int* channel,
    int* generator,
    int* switch_
);

void lusifer_density(
    double* g,
    int* generator,
    int* switch_
);

void lusifer_extra_max(
	int* maxex,
	int* maxgen
);

void lusifer_extra_data(
	int* gen,
	int* nch
);

void lusifer_extra_set(
	int* gen,
	int* nex,
	double* mw,
	double* gw,
	double* mz,
	double* gz,
	double* mh,
	double* gh,
	double* mt,
	double* gt
);

}
