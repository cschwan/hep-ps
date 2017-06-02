#ifdef HAVE_CONFIG_H
#include "config.hpp"
#endif

#define lusifer_initphasespace F77_FUNC_(lusifer_initphasespace, LUSIFER_INITPHASESPACE)
#define lusifer_phasespace F77_FUNC_(lusifer_phasespace, LUSIFER_PHASESPACE)
#define lusifer_density F77_FUNC_(lusifer_density, LUSIFER_DENSITY)

#define lusifer_extra_max F77_FUNC_(lusifer_extra_max, LUSIFER_EXTRA_MAX)
#define lusifer_extra_data F77_FUNC_(lusifer_extra_data, LUSIFER_EXTRA_DATA)
#define lusifer_extra_set F77_FUNC_(lusifer_extra_set, LUSIFER_EXTRA_SET)

constexpr std::size_t maxg = 1;
constexpr std::size_t maxe = 9;
constexpr std::size_t maxch = 20000;
constexpr std::size_t maxv = 40;

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

extern struct
{
	double alphaisr;
	double scale;
	double meisr;
	double s[1 << maxe];
	double p[1 << maxe][4];
	double mass[maxv + 1];
	double width[maxv + 1];
	int nchannel[maxg];
	int nexternal[maxg];
	int allbinary[maxg];
} lusifer_general_;

extern struct
{
	double powerinv[maxg][maxch][maxe];
	double mcutinv[maxg][(1 << maxe) + 1];
	int ininv[maxg][maxch][maxe];
	int idhepinv[maxg][maxch][maxv];
	int ninv[maxg][maxch];
	int lmin[maxg][maxch][maxe][maxe];
	int lmax[maxg][maxch][maxe][maxe];
} lusifer_cinv_;

extern struct
{
	int indecay[maxg][maxch][maxe];
	int out1decay[maxg][maxch][maxe];
	int out2decay[maxg][maxch][maxe];
	int ndecay[maxg][maxch];
} lusifer_cdecay_;

extern struct
{
	double powerprocess[maxg][maxch][maxe];
	double ccutprocess[maxg][1 << maxe];
	int in1process[maxg][maxch][maxe];
	int in2process[maxg][maxch][maxe];
	int out1process[maxg][maxch][maxe];
	int out2process[maxg][maxch][maxe];
	int inprocess[maxg][maxch][maxe];
	int virtprocess[maxg][maxch][maxe];
	int idhepprocess[maxg][maxch][maxe];
	int nprocess[maxg][maxch];
} lusifer_cprocess_;

extern struct
{
	int nsinv[maxg][maxch];
	int chinv[maxg][maxch];
	int maxinv[maxg];
	int ntprocess[maxg][maxch];
	int chprocess[maxg][maxch];
	int maxprocess[maxg];
	int nsdecay[maxg][maxch];
	int chdecay[maxg][maxch];
	int maxdecay[maxg];
	int numinv[maxg][maxch][maxe];
	int numprocess[maxg][maxch][maxe];
	int numdecay[maxg][maxch][maxe];
} lusifer_cdensity_;

}
