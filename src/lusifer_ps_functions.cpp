#include "hep/ps/lusifer_ps_functions.hpp"

#include <cmath>

namespace hep
{

template <typename T>
T lusifer_invariant_jacobian(T power, T mass, T width, T x, T xmin, T xmax)
{
    using std::atan;
    using std::copysign;
    using std::fabs;
    using std::log;
    using std::pow;

    T m2 = copysign(mass * mass, mass);
    mass = fabs(mass);

    if (width > T())
    {
        T const mg = mass * width;
        T const min = (xmin - m2) / mg;
        T const max = (xmax - m2) / mg;

        // formula for the difference of two arctans
        T max_min = atan((xmax - xmin) / mg / (T(1.0) + max * min));
        if (min * max < T(-1.0))
        {
            max_min += copysign(acos(T(-1.0)), max);
        }

        T const xprime = x - m2;

        return mg / (max_min * (xprime * xprime + mg * mg));
    }

    if (power == T())
    {
        return T(1.0) / (xmax - xmin);
    }

    if (mass == T())
    {
        m2 -= T(1e-6);
    }

    if (power == T(1.0))
    {
        // TODO: WARNING this branch is untested!
        return T(1.0) / (log((xmax - m2) / (xmin - m2)) * (x - m2));
    }

    T const omp = T(1.0) - power;

    return omp / ((pow(fabs(xmax - m2), omp) - pow(fabs(xmin - m2), omp)) *
        pow(fabs(x - m2), power));
}

template <typename T>
T lusifer_invariant_map(T power, T mass, T width, T x, T xmin, T xmax)
{
    using std::atan;
    using std::copysign;
    using std::exp;
    using std::fabs;
    using std::log;
    using std::pow;
    using std::tan;

    T m2 = copysign(mass * mass, mass);
    mass = fabs(mass);
    T result;

    if (width > T())
    {
        T const mg = mass * width;
        T const min = (xmin - m2) / mg;
        T const max = (xmax - m2) / mg;

        // formula for the difference of two arctans
        T max_min = atan((xmax - xmin) / mg / (T(1.0) + max * min));

        if (min * max < T(-1.0))
        {
            max_min += (max > T()) ? acos(T(-1.0)) : -acos(T(-1.0));
        }

        result = m2 + mg * tan(x * max_min + atan(min));
    }
    else if (power == T())
    {
        result = x * xmax + (T(1.0) - x) * xmin;

        // avoid infinite recursion
        if ((result < xmin) || (result > xmax))
        {
            // returning the upper bound produces the least problems
            return xmax;
        }
    }
    else
    {
        if (mass == T())
        {
            // this reduces the number of extremely small invariants which introduce numerical
            // problems in the matrix elements
            m2 -= T(1e-6);
        }

        if (power == T(1.0))
        {
            // TODO: WARNING this branch is untested!
            result = m2 + exp(x * log(xmax - m2) + (T(1.0) - x) * log(xmin - m2));
        }
        else
        {
            T const omp = T(1.0) - power;

            result = m2 + pow(x * pow(fabs(xmax - m2), omp) + (T(1.0) - x) *
                pow(fabs(xmin - m2), omp), T(1.0) / omp);
        }
    }

    // if the result is outside of the expected interval, calculate it linearly; this will only
    // happen if `xmax` is ridiculously small or if `xmin` is close to `xmax`, in which case it is a
    // valid approximation
    if ((result < xmin) || (result > xmax))
    {
        return lusifer_invariant_map(T(), T(), T(), x, xmin, xmax);
    }

    return result;
}

// -------------------- EXPLICIT TEMPLATE INSTANTIATIONS --------------------

template double lusifer_invariant_jacobian(double, double, double, double, double, double);

template double lusifer_invariant_map(double, double, double, double, double, double);

}
