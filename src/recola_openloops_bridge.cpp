#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <recola.h>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace
{

std::string const LO = "LO";
std::string const NLO = "NLO";

std::unordered_map<std::string, int> recola_process_list;
std::vector<int> recola_alphas_powers;
std::vector<int> recola_process_sizes;
std::vector<double> last_phase_space_point;

double muren = -1.0;
double alphas = -1.0;

template <typename T>
void ignore(T const&)
{
}

}

extern "C"
{

int order_ew = -1;
int loop_order_qcd = -1;

void ol_setparameter_int(const char* param, int val)
{
    if (std::strcmp(param, "ew_scheme") == 0)
    {
        // set alpha from Gmu, MW and MZ
        assert(val == 1);
    }
    else if (std::strcmp(param, "ew_renorm_scheme") == 0)
    {
        assert(val == 1);
    }
    else if (std::strcmp(param, "ew_renorm") == 0)
    {
    }
    else if (std::strcmp(param, "complex_mass_scheme") == 0)
    {
        assert(val == 2);

        Recola::set_complex_mass_scheme_rcl();
    }
    else if (std::strcmp(param, "polenorm") == 0)
    {
        if (val == 0)
        {
            // TODO: NYI
            assert( false );
        }
        else if (val == 1)
        {
            // here there's nothing to do
        }
        else
        {
            assert( false );
        }
    }
    else if (std::strcmp(param, "preset") == 0)
    {
    }
    else if (std::strcmp(param, "max_point") == 0)
    {
    }
    else if (std::strcmp(param, "redlib1") == 0)
    {
    }
    else if (std::strcmp(param, "redlib2") == 0)
    {
    }
    else if (std::strcmp(param, "use_coli_cache") == 0)
    {
    }
    else if (std::strcmp(param, "loop_order_ew") == 0)
    {
        // here we have nothing to do
    }
    else if (std::strcmp(param, "loop_order_qcd") == 0)
    {
        loop_order_qcd = val;
    }
    else if (std::strcmp(param, "order_ew") == 0)
    {
        order_ew = val;
    }
    else if (std::strcmp(param, "verbose") == 0)
    {
        Recola::set_print_level_squared_amplitude_rcl(val);
        Recola::set_print_level_correlations_rcl(val);
    }
    else if (std::strcmp(param, "write_parameters") == 0)
    {
    }
    else
    {
        std::cerr << param << '\n';
        assert(false);
    }
}

void ol_setparameter_double(const char* param, double val)
{
    static double mz = 0.0;
    static double gz = 0.0;
    static double mw = 0.0;
    static double gw = 0.0;
    static double mh = 0.0;
    static double gh = 0.0;
    static double mt = 0.0;
    static double gt = 0.0;

    if (std::strcmp(param, "mass(23)") == 0)
    {
        mz = val;

        if (mz != 0 && gz != 0)
        {
            Recola::set_pole_mass_z_rcl(mz, gz);
        }
    }
    else if (std::strcmp(param, "width(23)") == 0)
    {
        gz = val;

        if (mz != 0 && gz != 0)
        {
            Recola::set_pole_mass_z_rcl(mz, gz);
        }
    }
    else if (std::strcmp(param, "mass(24)") == 0)
    {
        mw = val;

        if (mw != 0 && gw != 0)
        {
            Recola::set_pole_mass_w_rcl(mw, gw);
        }
    }
    else if (std::strcmp(param, "width(24)") == 0)
    {
        gw = val;

        if (mw != 0 && gw != 0)
        {
            Recola::set_pole_mass_w_rcl(mw, gw);
        }
    }
    else if (std::strcmp(param, "mass(25)") == 0)
    {
        mh = val;

        if (mh != 0 && gh != 0)
        {
            Recola::set_pole_mass_h_rcl(mh, gh);
        }
    }
    else if (std::strcmp(param, "width(25)") == 0)
    {
        gh = val;

        if (mh != 0 && gh != 0)
        {
            Recola::set_pole_mass_h_rcl(mh, gh);
        }
    }
    else if (std::strcmp(param, "mass(6)") == 0)
    {
        mt = val;

        if (mt != 0 && gt != 0)
        {
            Recola::set_pole_mass_top_rcl(mt, gt);
        }
    }
    else if (std::strcmp(param, "width(6)") == 0)
    {
        gt = val;

        if (mt != 0 && gt != 0)
        {
            Recola::set_pole_mass_top_rcl(mt, gt);
        }
    }
    else if (std::strcmp(param, "mass(5)") == 0)
    {
    }
    else if (std::strcmp(param, "width(5)") == 0)
    {
    }
    else if (std::strcmp(param, "mass(15)") == 0)
    {
    }
    else if (std::strcmp(param, "gmu") == 0)
    {
        Recola::use_gfermi_scheme_real_and_set_gfermi_rcl(val);
    }
    else if (std::strcmp(param, "alphas") == 0)
    {
        alphas = val;

        if (muren != -1.0)
        {
            Recola::set_alphas_rcl(alphas, muren, 5);
        }
        else
        {
            Recola::set_alphas_rcl(alphas, -1.0, 5);
        }
    }
    else if (std::strcmp(param, "mureg") == 0)
    {
        Recola::set_mu_ir_rcl(val);
        Recola::set_mu_uv_rcl(val);
    }
    else if (std::strcmp(param, "muren") == 0)
    {
        muren = val;

        if (alphas != -1.0 && muren != -1.0)
        {
            Recola::set_alphas_rcl(alphas, muren, 5);
            alphas = -1.0;
            muren = -1.0;
        }
    }
    else
    {
        std::cerr << param << '\n';
        assert(false);
    }
}

void ol_setparameter_string(const char* param, char* val)
{
    ignore(param);
    ignore(val);

    assert(false);
}

void ol_getparameter_int(const char* param, int* val)
{
    ignore(param);
    ignore(val);

    assert(false);
}

void ol_getparameter_double(const char* param, double* val)
{
    if (std::strcmp(param, "alphas") == 0)
    {
        Recola::get_alphas_rcl(*val);
    }
    else
    {
        std::cerr << param << '\n';
        assert(false);
    }
}

int ol_register_process(const char* process, int amptype)
{
    std::string recola_amptype;

    switch (amptype)
    {
    case 1:
        recola_amptype = "LO";
        break;

    case 11:
        recola_amptype = "NLO";
        break;

    default:
        std::cerr << amptype << '\n';
        assert(false);
    }

    int final_states = -2;

    std::istringstream input(process);
    std::string recola_process;

    for (std::string token; std::getline(input, token, ' '); )
    {
        if (token == "->")
        {
            recola_process.append(token);
            recola_process.append(" ");

            continue;
        }

        switch (std::stoi(token))
        {
        case -13:
           recola_process.append("mu+ ");
           break;

        case -11:
           recola_process.append("e+ ");
           break;

        case -4:
            recola_process.append("c~ ");
            break;

        case -3:
            recola_process.append("s~ ");
            break;

        case -2:
            recola_process.append("u~ ");
            break;

        case -1:
            recola_process.append("d~ ");
            break;

        case 1:
            recola_process.append("d ");
            break;

        case 2:
            recola_process.append("u ");
            break;

        case 3:
            recola_process.append("s ");
            break;

        case 4:
            recola_process.append("c ");
            break;

        case 12:
           recola_process.append("nu_e ");
           break;

        case 14:
           recola_process.append("nu_mu ");
           break;

        case 21:
           recola_process.append("g ");
           break;

        default:
            std::cerr << token << '\n';
            assert(false);
        }

        ++final_states;
    }

    auto const result = recola_process_list.find(recola_process);

    if (result == recola_process_list.end())
    {
        int process_id = recola_process_list.size() + 1;

        Recola::define_process_rcl(process_id, recola_process.c_str(), recola_amptype.c_str());

        int as_power = -1;

        if (recola_amptype == NLO)
        {
            Recola::unselect_all_gs_powers_LoopAmpl_rcl(process_id);
            Recola::unselect_all_gs_powers_BornAmpl_rcl(process_id);

            // Recola::select_gs_power_BornAmpl_rcl(1, alphas);

            // TODO: this is process-specific information
            std::array<int, 2> born_gs = { 0, 2 };
            std::array<int, 3> loop_gs = { 0, 2, 4 };

            as_power = loop_order_qcd;

            for (int born : born_gs)
            {
                for (int loop : loop_gs)
                {
                    if ((born + loop) == 2 * loop_order_qcd)
                    {
                        Recola::select_gs_power_LoopAmpl_rcl(process_id, loop);
                        Recola::select_gs_power_BornAmpl_rcl(process_id, born);
                    }
                }
            }
        }
        else
        {
            Recola::unselect_all_gs_powers_BornAmpl_rcl(process_id);

            // TODO: this is process-specific information
            std::array<int, 2> born_gs = { 0, 2 };

            as_power = final_states - order_ew;

            for (int born1 : born_gs)
            {
                for (int born2 : born_gs)
                {
                    if ((born1 + born2) == 2 * as_power)
                    {
                        Recola::select_gs_power_BornAmpl_rcl(process_id, born1);
                        Recola::select_gs_power_BornAmpl_rcl(process_id, born2);
                    }
                }
            }
        }

        recola_process_list.emplace(recola_process, process_id);
        recola_alphas_powers.push_back(as_power);
        recola_process_sizes.push_back(final_states + 2);

        return process_id++;
    }
    else
    {
        std::cerr << "found process: " << result->second << '\n';

        return result->second;
    }
}

int ol_n_external(int id)
{
    ignore(id);

    assert(false);
}

void ol_start()
{
    Recola::generate_processes_rcl();
}

void ol_finish()
{
}

void ol_evaluate_tree(int id, double* pp, double* m2tree)
{
    int n = recola_process_sizes.at(id - 1);

    last_phase_space_point.resize(4 * n);

    for (int i = 0; i != n; ++i)
    {
        last_phase_space_point.at(4 * i + 0) = pp[5 * i + 0];
        last_phase_space_point.at(4 * i + 1) = pp[5 * i + 1];
        last_phase_space_point.at(4 * i + 2) = pp[5 * i + 2];
        last_phase_space_point.at(4 * i + 3) = pp[5 * i + 3];
        // ignore mass
    }

    auto mom = reinterpret_cast <double (*) [4]> (last_phase_space_point.data());

    Recola::compute_process_rcl(id, mom, LO);
    int alphas = recola_alphas_powers.at(id - 1);
    Recola::get_squared_amplitude_rcl(id, alphas, LO, *m2tree);
}

void ol_evaluate_cc(int id, double* pp, double* m2tree, double* m2cc, double* m2ew)
{
    *m2tree = 0.0;
    *m2ew = 0.0;

    int n = recola_process_sizes.at(id - 1);

    last_phase_space_point.resize(4 * n);

    for (int i = 0; i != n; ++i)
    {
        last_phase_space_point.at(4 * i + 0) = pp[5 * i + 0];
        last_phase_space_point.at(4 * i + 1) = pp[5 * i + 1];
        last_phase_space_point.at(4 * i + 2) = pp[5 * i + 2];
        last_phase_space_point.at(4 * i + 3) = pp[5 * i + 3];
        // ignore mass
    }

    auto mom = reinterpret_cast <double (*) [4]> (last_phase_space_point.data());

    for (int i = 0; (i + 1) < n; ++i)
    {
        for (int j = i + 1; j != n; ++j)
        {
            double& a2cc = m2cc[i + j * (j - 1) / 2];
            Recola::compute_colour_correlation_rcl(id, mom, i + 1, j + 1);
            // TODO: generalize alphas=1
            Recola::get_colour_correlation_rcl(id, 1, i + 1, j + 1, a2cc);
            // TODO: this is dependent on the particle type
            a2cc *= 4.0 / 3.0;
        }
    }

    Recola::compute_process_rcl(id, mom, LO);
    // TODO: generalize alphas=1
    Recola::get_squared_amplitude_rcl(id, 1, LO, *m2tree);
}

void ol_evaluate_sc(int id, double* pp, int emitter, double* polvect, double* m2sc)
{
    ignore(id);
    ignore(pp);
    ignore(emitter);
    ignore(polvect);
    ignore(m2sc);

    assert(false);
}

void ol_evaluate_loop(int id, double* pp, double* m2tree, double* m2loop, double* acc)
{
    int n = recola_process_sizes.at(id - 1);

    last_phase_space_point.resize(4 * n);

    for (int i = 0; i != n; ++i)
    {
        last_phase_space_point.at(4 * i + 0) = pp[5 * i + 0];
        last_phase_space_point.at(4 * i + 1) = pp[5 * i + 1];
        last_phase_space_point.at(4 * i + 2) = pp[5 * i + 2];
        last_phase_space_point.at(4 * i + 3) = pp[5 * i + 3];
        // ignore mass
    }

    auto mom = reinterpret_cast <double (*) [4]> (last_phase_space_point.data());

    double result[2] = {};
    Recola::compute_process_rcl(id, mom, "NLO", result);
    int alphas = recola_alphas_powers.at(id - 1);
    Recola::get_squared_amplitude_rcl(id, alphas, "NLO", result[1]);

    m2loop[0] = result[1];
    m2loop[1] = -1.0;
    m2loop[2] = -1.0;

    // TODO: this isn't the tree-level result, but we need to escape the k-factor tech cut
    *m2tree = result[1];
}

void ol_evaluate_loop2(int id, double* pp, double* m2loop2, double* acc)
{
    ignore(id);
    ignore(pp);
    ignore(m2loop2);
    ignore(acc);

    assert(false);
}

void ol_evaluate_ct(int id, double* pp, double* m2_tree, double* m2_ct)
{
    ignore(id);
    ignore(pp);
    ignore(m2_tree);
    ignore(m2_ct);

    // TODO: for the time being we don't need scale variations
    *m2_ct = 0.0;
}

void ol_evaluate_full(int id, double* pp, double* m2tree, double* m2loop, double* m2ir, double* m2loop2, double* m2iop, double* acc)
{
    ignore(id);
    ignore(pp);
    ignore(m2tree);
    ignore(m2loop);
    ignore(m2ir);
    ignore(m2loop2);
    ignore(m2iop);
    ignore(acc);

    assert(false);
}

}
