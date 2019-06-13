/**
 * @file harmonic_oscillator.hpp
 * @brief Simple `HarmonicOscillator` class is defined for the code test.
 */
#ifndef _HARMONIC_OSCILLATOR_H_
#define _HARMONIC_OSCILLATOR_H_

#include "model.hpp"



/**
 * @class HarmonicOscillator
 * @brief 
 */
class HarmonicOscillator final: public Model
{
    public:
        HarmonicOscillator     ( std::string name ): Model(name) {}
        /*
        double V   ( double** f, int n, int idx, double a=1 ) const override
        { return pow(f[n][idx],2) / (2*a*a); }
        double dV  ( double** f, int n, int idx, double a=1 ) const override
        { return f[n][idx]/a; }
        */
};



#endif