#ifndef _CREATE_EVOLUTION_SCHEME_H_
#define _CREATE_EVOLUTION_SCHEME_H_

#include "leapfrog_in_Minkowsky.hpp"
#include "leapfrog_in_FRW.hpp"



template <typename T=EvolutionScheme>
std::unique_ptr<T> createEvolutionScheme (std::string name, int precision, Expansion expansion,
                                          double** f, double** df)
{
    if ( name == "leapfrog with ABC" ) return std::make_unique<LeapFrogWithABC>(precision);
    if ( name == "leapfrog" ){
        if( expansion == Expansion::no_expansion )
            return std::make_unique<LeapFrogInMinkowsky>(precision);
        else if( expansion == Expansion::self_consist )
            return std::make_unique<LeapFrogSelf>(precision, f, df);
        else
            return std::make_unique<LeapFrogRadMat>(precision);
    }
    return nullptr;
}


#endif