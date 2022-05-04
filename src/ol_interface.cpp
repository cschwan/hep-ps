#include "hep/ps/ol_interface.hpp"
#include "hep/ps/suppress_banners.hpp"

#include "config.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
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

namespace
{

template <typename T>
void ignore(T const&)
{
}

std::string normalize_string(char const* string)
{
    std::stringstream ss(string);
    std::string result;
    std::string s;

    while (std::getline(ss, s, ' '))
    {
        //std::cout << "}}} " << s << '\n';

        if (!s.empty())
        {
            result.append(s);
            result.append(" ");
        }
    }

    //std::cout << "normalize " << string << " to " << result << '\n';

    return result;
}

}

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
    : started_(false)
    , set_order_qcd_{false}
    , replacement_rules_()
    , zero_rules_()
{
#ifdef HAVE_OPENLOOPS
    if (suppress_banners())
    {
        ol_setparameter_int("nosplash", 1);
    }
#endif
}

ol_interface::~ol_interface()
{
#ifdef HAVE_OPENLOOPS
    if (started_)
    {
        ol_finish();
    }
#endif
}

void ol_interface::setparameter_int(char const* param, int val)
{
#ifdef HAVE_OPENLOOPS
    ol_setparameter_int(param, val);
#else
    ignore(param);
    ignore(val);
#endif
}

void ol_interface::setparameter_double(char const* param, double val)
{
#ifdef HAVE_OPENLOOPS
    ol_setparameter_double(param, val);
#else
    ignore(param);
    ignore(val);
#endif
}

void ol_interface::setparameter_string(char const* param, char* val)
{
#ifdef HAVE_OPENLOOPS
    ol_setparameter_string(param, val);
#else
    ignore(param);
    ignore(val);
#endif
}

void ol_interface::getparameter_int(char const* param, int* val)
{
#ifdef HAVE_OPENLOOPS
    ol_getparameter_int(param, val);
#else
    ignore(param);
    ignore(val);
#endif
}

void ol_interface::getparameter_double(char const* param, double* val)
{
#ifdef HAVE_OPENLOOPS
    ol_getparameter_double(param, val);
#else
    ignore(param);
    ignore(val);
#endif
}

int ol_interface::register_process(char const* process, int amptype)
{
#ifdef HAVE_OPENLOOPS
    return ol_register_process(process, amptype);
#else
    ignore(process);
    ignore(amptype);

    return -1;
#endif
}

int ol_interface::n_external(int id)
{
#ifdef HAVE_OPENLOOPS
    return ol_n_external(id);
#else
    ignore(id);

    throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_tree(int id, double* pp, double* m2tree)
{
#ifdef HAVE_OPENLOOPS
    if (!started_)
    {
        started_ = true;
        ol_start();
    }

    if (id == 0)
    {
        *m2tree = 0.0;
    }
    else
    {
        ol_evaluate_tree(id, pp, m2tree);
    }
#else
    ignore(started_);
    ignore(id);
    ignore(pp);
    ignore(m2tree);

    throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_cc(int id, double* pp, double* m2tree, double* m2cc, double* m2ew)
{
#ifdef HAVE_OPENLOOPS
    if (!started_)
    {
        started_ = true;
        ol_start();
    }

    if (id == 0)
    {
        *m2tree = 0.0;
        *m2cc = 0.0;
        *m2ew = 0.0;
    }
    else
    {
        ol_evaluate_cc(id, pp, m2tree, m2cc, m2ew);
    }
#else
    ignore(started_);
    ignore(id);
    ignore(pp);
    ignore(m2tree);
    ignore(m2cc);
    ignore(m2ew);

    throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_sc(int id, double* pp, int emitter, double* polvect, double* m2sc)
{
#ifdef HAVE_OPENLOOPS
    if (!started_)
    {
        started_ = true;
        ol_start();
    }

    if (id == 0)
    {
        *m2sc = 0.0;
    }
    else
    {
        ol_evaluate_sc(id, pp, emitter, polvect, m2sc);
    }
#else
    ignore(started_);
    ignore(id);
    ignore(pp);
    ignore(emitter);
    ignore(polvect);
    ignore(m2sc);

    throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_loop(int id, double* pp, double* m2tree, double* m2loop, double* acc)
{
#ifdef HAVE_OPENLOOPS
    if (!started_)
    {
        started_ = true;
        ol_start();
    }

    if (id == 0)
    {
        *m2tree = 0.0;
        *m2loop = 0.0;
        *acc = 0.0;
    }
    else
    {
        ol_evaluate_loop(id, pp, m2tree, m2loop, acc);
    }
#else
    ignore(started_);
    ignore(id);
    ignore(pp);
    ignore(m2tree);
    ignore(m2loop);
    ignore(acc);

    throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_loop2(int id, double* pp, double* m2loop2, double* acc)
{
#ifdef HAVE_OPENLOOPS
    if (!started_)
    {
        started_ = true;
        ol_start();
    }

    if (id == 0)
    {
        *m2loop2 = 0.0;
        *acc = 0.0;
    }
    else
    {
        ol_evaluate_loop2(id, pp, m2loop2, acc);
    }
#else
    ignore(started_);
    ignore(id);
    ignore(pp);
    ignore(m2loop2);
    ignore(acc);

    throw std::runtime_error("OpenLoops support not enabled");
#endif
}

void ol_interface::evaluate_ct(int id, double* pp, double* m2_tree, double* m2_ct)
{
#ifdef HAVE_OPENLOOPS
    if (!started_)
    {
        started_ = true;
        ol_start();
    }

    if (id == 0)
    {
        *m2_tree = 0.0;
        *m2_ct = 0.0;
    }
    else
    {
        ol_evaluate_ct(id, pp, m2_tree, m2_ct);
    }
#else
    ignore(started_);
    ignore(id);
    ignore(pp);
    ignore(m2_tree);
    ignore(m2_ct);

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
    if (!started_)
    {
        started_ = true;
        ol_start();
    }

    if (id == 0)
    {
        *m2tree = 0.0;
        *m2loop = 0.0;
        *m2ir1 = 0.0;
        *m2loop2 = 0.0;
        *m2ir2 = 0.0;
        *acc = 0.0;
    }
    else
    {
        ol_evaluate_full(id, pp, m2tree, m2loop, m2ir1, m2loop2, m2ir2, acc);
    }
#else
    ignore(started_);
    ignore(id);
    ignore(pp);
    ignore(m2tree);
    ignore(m2loop);
    ignore(m2ir1);
    ignore(m2loop2);
    ignore(m2ir2);
    ignore(acc);

    throw std::runtime_error("OpenLoops support not enabled");
#endif
}

int ol_interface::register_process(char const* process, me_type type, int order_qcd, int order_ew)
{
    std::string p = normalize_string(process);

    if (std::any_of(zero_rules_.begin(), zero_rules_.end(),
        [=](std::tuple<std::string, int, int> const& rule) {
            return (p == std::get<0>(rule)) && (order_qcd == std::get<1>(rule)) &&
                (order_ew == std::get<2>(rule));
        }))
    {
        //std::cout << ">>> replaced process: " << process << " at order O(as^" << order_qcd
        //    << ") with zero\n";

        return 0;
    }

    auto const& tuple = std::make_tuple(p, order_qcd, order_ew);
    auto const& rule = replacement_rules_.find(tuple);

    if (rule != replacement_rules_.end())
    {
        p = rule->second;
        //std::cout << ">>> replaced process: " << process << " with " << p << " at order O(as^"
        //    << order_qcd << ")\n";
    }
    else
    {
        //std::cout << ">>> didn't replace process: " << process << " at order O(as^" << order_qcd
        //    << ")\n";
    }

    int amptype = 1;

    switch (type)
    {
    case me_type::born:
        setparameter_int("loop_order_ew", -1);
        setparameter_int("loop_order_qcd", -1);

        if (set_order_qcd_)
        {
            setparameter_int("order_qcd", order_qcd);
        }
        else
        {
            setparameter_int("order_ew", order_ew);
        }

        amptype = 1;

        break;

    case me_type::color_correlated:
    case me_type::spin_correlated:
    case me_type::loop:
        setparameter_int("loop_order_ew", order_ew);
        setparameter_int("loop_order_qcd", order_qcd);

        amptype = 11;

        break;

    default:
        assert( false );
    }

    int const result = register_process(p.c_str(), amptype);

    // we'll reserve ID=0 for zero matrix elements
    assert( result != 0 );

    if (result == -1)
    {
        std::cerr << "couldn't find process: " << p << '\n';
        std::cerr << "amptype: " << amptype << '\n';
        std::cerr << "qcd: " << order_qcd << '\n';
        std::cerr << "ew: " << order_ew << '\n';
        //std::exit(1);
    }

    return result;
}

void ol_interface::register_replacement_rule(
    char const* process,
    char const* replacement,
    int order_qcd,
    int order_ew
) {
    replacement_rules_.emplace(std::make_tuple(normalize_string(process), order_qcd, order_ew),
        normalize_string(replacement));
}

void ol_interface::register_zero_rule(
    char const* process,
    int order_qcd,
    int order_ew
) {
    zero_rules_.emplace_back(std::make_tuple(normalize_string(process), order_qcd, order_ew));
}

}
