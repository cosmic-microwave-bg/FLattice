/**
 * @file model.hpp
 * @brief class `Model` is written where the potential is defined.
 */
#ifndef _MODEL_H_
#define _MODEL_H_

#include <cmath>
#include "parameter.hpp"



/**
 * @class Model
 */
class Model
{
    protected:
        std::string _name;

    public:
        Model      ( std::string name ): _name(name) {}
        
        /**
         * @brief   V, dV
         * @details When you simulate in the Minkowski spacetime, just write the potential \f$ V(\phi) \f$.
         *          When you simulate in the FRW spacetime, write the rescaled potential \f$ V(\bar{\phi}) \f$
         *          where \f$ \bar{\phi} = \phi / a \f$
         */
        virtual double V    ( double** f, int n, int idx, double a=1 ) const = 0;
        virtual double dV   ( double** f, int n, int idx, double a=1 ) const = 0;
        std::string    name ()                                         const { return _name; }
};



#endif