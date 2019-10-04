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
 * @details When you simulate in the Minkowski spacetime, just write the potential \f$ V(\phi) \f$.
 *          When you simulate in the FRW spacetime, write the rescaled potential \f$ V(\bar{\phi}) \f$
 *          where \f$ \bar{\phi} = \phi / a \f$
 */
template <typename T = double>
T V    (T** f, int n, int idx)
{ return pow(f[n][idx], 2) / 2; }

template <typename T = double>
T dV   (T** f, int n, int idx)
{ return f[n][idx]; }


#endif