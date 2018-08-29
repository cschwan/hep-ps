#ifdef HAVE_CONFIG_H
#include "config.hpp"
#endif

#define cofferaa_initgenerator F77_FUNC_(cofferaa_initgenerator, COFFERAA_INITGENERATOR)
#define cofferaa_generation F77_FUNC_(cofferaa_generation, COFFERAA_GENERATION)
#define cofferaa_density F77_FUNC_(cofferaa_density, COFFERAA_DENSITY)

#define cofferaa_extra_max F77_FUNC_(cofferaa_extra_max, COFFERAA_EXTRA_MAX)
#define cofferaa_extra_data F77_FUNC_(cofferaa_extra_data, COFFERAA_EXTRA_DATA)
#define cofferaa_extra_set F77_FUNC_(cofferaa_extra_set, COFFERAA_EXTRA_SET)

extern "C"
{

void cofferaa_generation(
    double const* random,
    double* kbeam,
    double* k,
    double* g,
    int* channel,
    int* generator,
    int* switch_
);

void cofferaa_density(
    double* k,
    double* g,
    int* generator,
    int* switch_
);

void cofferaa_initgenerator(
    double* energy,
    double* smin,
    int const* hepnum,
    int* generator,
    int* next,
    int* smodel,
    int* sincludecuts,
    int* ssub,
    int*,
    int*,
    int*,
    int*
);

void cofferaa_extra_max(
    int* maxex,
    int* maxgen
);

void cofferaa_extra_data(
    int* gen,
    int* nch
);

void cofferaa_extra_set(
    int* gen,
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
