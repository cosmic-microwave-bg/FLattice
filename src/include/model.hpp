/**
 * @file model.hpp
 * @brief `Model` is written where the potential is defined.
 */
#ifndef _MODEL_H_
#define _MODEL_H_

#include <cmath>
#include "parameter.hpp"


/**
 * @brief   V, dV
 * @details When you simulate in either spacetime, just write down the potential \f$ V(\phi) \f$.
 */
template <typename T = double>
T V    (T** f, int n, int idx)
{ return pow(f[n][idx], 2) / 2; }

template <typename T = double>
T dV   (T** f, int n, int idx)
{ return f[n][idx]; }


#endif