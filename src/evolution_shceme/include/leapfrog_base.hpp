/**
 * @file leapfrog_base.hpp
 * @brief The base class for the leapfrog evoltuion derived by EvolutioinScheme.
 */
#ifndef _LEAPFROG_BASE_H_
#define _LEAPFROG_BASE_H_

#include "evolution_base.hpp"



class LeapFrogBase: public EvolutionScheme
{
    protected:
        int _precision;
        /** @brief           The function to calculate the laplacian by fourth-order centeral differencial scheme.
         *  @param[in]  *f   The argument of the fields must be the pointer type for the vector optimization of the intel compiler!
         */
        double   laplacian        ( const double* f, int i, int j=0, int k=0 ) const;
        void     evolFields       ( double** f, double** df, const double h );
        void     evolFieldDerivs  ( double** f, double** df, const double h, double a=1 );

    public:
        LeapFrogBase              ( std::shared_ptr<Model> model, int precision )
        : EvolutionScheme(model), _precision(precision) {}
        virtual  ~LeapFrogBase    () {};
};



#endif