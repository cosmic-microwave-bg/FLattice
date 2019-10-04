/**
 * @file evolution_base.hpp
 * @brief The base class necessary for the time evolution is defined.
 */
#ifndef _EVOULUTION_SCHEME_H_
#define _EVOULUTION_SCHEME_H_

#include <memory>
#include "model.hpp"


/** 
 * @class EvolutionScheme
 * @brief The base class of time evolution of the fields.
 */
class EvolutionScheme
{
    protected:
        int _precision;

    public:
        EvolutionScheme          (int precision): _precision(precision) {}

        virtual ~EvolutionScheme () {}
        virtual void evolution   (double** f, double** df, int ouptut_step) = 0;
};



#endif