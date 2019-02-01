#ifndef HEP_PS_OL_INTERFACE_HPP
#define HEP_PS_OL_INTERFACE_HPP

/*
 * hep-ps - A C++ Library of Phase Space Integrands for High Energy Physics
 * Copyright (C) 2018-2019  Christopher Schwan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hep/ps/me_type.hpp"

namespace hep
{

/// Singleton for accessing the OpenLoops Matrix elements. If support for OpenLoops in this library
/// wasn't activated, some of the methods frow an exception. Currently this interface is only a thin
/// wrapper which automatically calls `ol_start` and `ol_finish`.
class ol_interface final
{
public:
    /// Returns the singleton instance.
    static ol_interface& instance();

    /// Returns `true` if this library was compiled against OpenLoops. Otherwise this function
    /// returns `false` and the member functions of this class won't work.
    static bool enabled();

    /// Destructor.
    ~ol_interface();

    /// Set integer parameter.
    void setparameter_int(char const* param, int val);

    /// Set double parameter.
    void setparameter_double(char const* param, double val);

    /// Set string parameter.
    void setparameter_string(char const* param, char* val);

    /// Return the value of an integer parameter.
    void getparameter_int(char const* param, int* val);

    /// Return the value of a double parameter.
    void getparameter_double(char const* param, double* val);

    /// Register a partonic process.
    int register_process(char const* process, int amptype);

    /// Query the number of external particle for a given process.
    int n_external(int id);

    /// Evaluate a tree matrix element.
    void evaluate_tree(int id, double* pp, double* m2tree);

    /// Evaluate a color-correlated matrix element.
    void evaluate_cc(int id, double* pp, double* m2tree, double* m2cc, double* m2ew);

    /// Evaluate a spin-correlated matrix element.
    void evaluate_sc(int id, double* pp, int emitter, double* polvect, double* m2sc);

    /// Evaluate a loop matrix element.
    void evaluate_loop(int id, double* pp, double* m2tree, double* m2loop, double* acc);

    /// Evaluate a squared loop matrix element.
    void evaluate_loop2(int id, double* pp, double* m2loop2, double* acc);

    /// Evaluate counterterms.
    void evaluate_ct(int id, double* pp, double* m2_tree, double* m2_ct);

    /// Evaluate everything.
    void evaluate_full(
        int id,
        double* pp,
        double* m2tree,
        double* m2loop,
        double* m2ir1,
        double* m2loop2,
        double* m2ir2,
        double* acc
    );

    /// Register a matrix element for the process represented by `process`, of the given `type` and
    /// coupling orders given by `order_qcd` and `order_ew`. The return value denotes an ID which is
    /// needed as an argument for other functions.
    int register_process(char const* process, me_type type, int order_qcd, int order_ew);

private:
    bool started_;
    bool set_order_qcd_;

    ol_interface();
};

}

#endif
