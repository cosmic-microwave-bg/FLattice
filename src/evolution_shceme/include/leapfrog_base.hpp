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
        /** @brief           The function to calculate the laplacian by fourth-order centeral differencial scheme.
         *  @param[in]  *f   The argument of the fields must be the pointer type for the vector optimization of the intel compiler!
         */
        const double _C[4] = {+0.675603595979828817023844, -0.175603595979828817023844, -0.175603595979828817023844, +0.675603595979828817023844};
        const double _D[3] = {+1.351207191959657634047688, -1.702414383919315268095376, +1.351207191959657634047688};

        double   laplacian        (const double* f, int i, int j=0, int k=0) const;
        void     evolFields       (double** f, double** df, const double h, double a=1);
        void     evolFieldDerivs  (double** f, double** df, const double h, double a=1);

    public:
        LeapFrogBase              (int precision)
        : EvolutionScheme(precision) {}
        
        virtual  ~LeapFrogBase    () {};
};



#endif