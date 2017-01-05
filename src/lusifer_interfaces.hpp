#ifdef HAVE_CONFIG_H
#include "config.hpp"
#endif

#define lusifer_initphasespace F77_FUNC(lusifer_initphasespace, INITPHASESPACE)
#define lusifer_phasespace F77_FUNC(lusifer_phasespace, PHASESPACE)
#define lusifer_density F77_FUNC(lusifer_density, DENSITY)

#define lusifer_extra_generatormax F77_FUNC(generatormax, GENERATORMAX)
#define lusifer_extra_generatordata F77_FUNC(generatordata, GENERATORDATA)
#define lusifer_extra_generatorset F77_FUNC(generatorset, GENERATORSET)

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

void lusifer_extra_generatormax(
	int* maxex,
	int* maxgen
);

void lusifer_extra_generatordata(
	int* gen,
	int* nch
);

void lusifer_extra_generatorset(
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
