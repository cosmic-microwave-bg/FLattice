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
{ return pow(f[n][idx], 2) / 2; }  // harmonic oscillation
// { return (1 - pow(1 + f[n][idx]*f[n][idx], -p)) / (2*p); }  // ULAP potential
// { return ( 1 + K*log((pow(f[0][idx],2) + pow(f[1][idx],2))/2) ) * pow(f[n][idx],2)/2; }  // gravity mediated Q-ball potential


template <typename T = double>
T dV   (T** f, int n, int idx)
{ return f[n][idx]; }  // harmonic oscillation
// { return pow(1 + f[n][idx]*f[n][idx], -p-1) * f[n][idx]; }  // ULAP potential
// { return ( 1+K + K*log((pow(f[0][idx],2) + pow(f[1][idx],2))/2) ) * f[n][idx]; }  // gravity mediated Q-ball potential


#endif