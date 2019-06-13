#ifndef _CREATE_MODEL_H_
#define _CREATE_MODEL_H_

#include "harmonic_oscillator.hpp"
#include "qball.hpp"



template <typename T=Model>
std::shared_ptr<T> createModel (std::string name)
{
    //if( name == "harmonic_oscillator" ) return std::make_shared<HarmonicOscillator>(name);
    if( name == "Qball"               ) return std::make_shared<Qball>(name);

    return nullptr;
}



#endif