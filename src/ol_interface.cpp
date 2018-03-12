#include "hep/ps/ol_interface.hpp"

#include "config.hpp"

#include <stdexcept>

#ifdef HAVE_OPENLOOPS
extern "C"
{

void ol_setparameter_int(const char* param, int val);
void ol_setparameter_double(const char* param, double val);
void ol_setparameter_string(const char* param, char* val);
void ol_getparameter_int(const char* param, int* val);
void ol_getparameter_double(const char* param, double* val);
int ol_register_process(const char* process, int amptype);
int ol_n_external(int id);
//void ol_phase_space_point(int id, double sqrt_s, double* pp);
void ol_start();
void ol_finish();
void ol_evaluate_tree(int id, double* pp, double* m2tree);
void ol_evaluate_cc(int id, double* pp, double* m2tree, double* m2cc, double* m2ew);
void ol_evaluate_sc(int id, double* pp, int emitter, double* polvect, double* m2sc);
void ol_evaluate_loop(int id, double* pp, double* m2tree, double* m2loop, double* acc);
void ol_evaluate_loop2(int id, double* pp, double* m2loop2, double* acc);
void ol_evaluate_ct(int id, double* pp, double* m2_tree, double* m2_ct);
void ol_evaluate_full(int id, double* pp, double* m2tree, double* m2loop, double* m2ir, double* m2loop2, double* m2iop, double* acc);

}
#endif

namespace hep
{

ol_interface& ol_interface::instance()
{
	static ol_interface object;
	return object;
}

bool ol_interface::enabled()
{
#ifdef HAVE_OPENLOOPS
	return true;
#else
	return false;
#endif
}

ol_interface::ol_interface()
{
#ifdef HAVE_OPENLOOPS
	ol_start();
#endif
}

ol_interface::~ol_interface()
{
#ifdef HAVE_OPENLOOPS
	ol_finish();
#endif
}

void ol_interface::setparameter_int(char const* param, int val)
{
#ifdef HAVE_OPENLOOPS
	ol_setparameter_int(param, val);
#endif
}

void ol_interface::setparameter_double(char const* param, double val)
{
#ifdef HAVE_OPENLOOPS
	ol_setparameter_double(param, val);
#endif
}

void ol_interface::setparameter_string(char const* param, char* val)
{
#ifdef HAVE_OPENLOOPS
	ol_setparameter_string(param, val);
#endif
}

void ol_interface::getparameter_int(char const* param, int* val)
{
#ifdef HAVE_OPENLOOPS
	ol_getparameter_int(param, val);
#endif
}

void ol_interface::getparameter_double(char const* param, double* val)
{
#ifdef HAVE_OPENLOOPS
	ol_getparameter_double(param, val);
#endif
}

int ol_interface::register_process(char const* process, int amptype)
{
#ifdef HAVE_OPENLOOPS
	ol_register_process(process, amptype);
#endif
}

int ol_interface::n_external(int id)
{
#ifdef HAVE_OPENLOOPS
	ol_n_external(id);
#else
	throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_tree(int id, double* pp, double* m2tree)
{
#ifdef HAVE_OPENLOOPS
	ol_evaluate_tree(id, pp, m2tree);
#else
	throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_cc(
	int id,
	double* pp,
	double* m2tree,
	double* m2cc,
	double* m2ew
) {
#ifdef HAVE_OPENLOOPS
	ol_evaluate_cc(id, pp, m2tree, m2cc, m2ew);
#else
	throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_sc(
	int id,
	double* pp,
	int emitter,
	double* polvect,
	double* m2sc
) {
#ifdef HAVE_OPENLOOPS
	ol_evaluate_sc(id, pp, emitter, polvect, m2sc);
#else
	throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_loop(
	int id,
	double* pp,
	double* m2tree,
	double* m2loop,
	double* acc
) {
#ifdef HAVE_OPENLOOPS
	ol_evaluate_loop(id, pp, m2tree, m2loop, acc);
#else
	throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_loop2(
	int id,
	double* pp,
	double* m2loop2,
	double* acc
) {
#ifdef HAVE_OPENLOOPS
	ol_evaluate_loop2(id, pp, m2loop2, acc);
#else
	throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_ct(
	int id,
	double* pp,
	double* m2_tree,
	double* m2_ct
) {
#ifdef HAVE_OPENLOOPS
	ol_evaluate_ct(id, pp, m2_tree, m2_ct);
#else
	throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_full(
	int id,
	double* pp,
	double* m2tree,
	double* m2loop,
	double* m2ir1,
	double* m2loop2,
	double* m2ir2,
	double* acc
) {
#ifdef HAVE_OPENLOOPS
	ol_evaluate_full(id, pp, m2tree, m2loop, m2ir1, m2loop2, m2ir2, acc);
#else
	throw std::runtime_error("OpenLoops support not enabled");
#endif
}

}
