#ifndef _LEAPFROG_WITH_EXPANSION_H_
#define _LEAPFROG_WITH_EXPANSION_H_

#include <cmath>
#include "leapfrog_base.hpp"
#include "calculation.hpp"



class LeapFrogRadMat final: public LeapFrogBase
{
    public:
        LeapFrogRadMat  ( std::shared_ptr<Model> model, int precision )
        : LeapFrogBase(model, precision) 
        {
            if( expansion == Expansion::rad_domination ) da = 1/t0;
            if( expansion == Expansion::mat_domination ) da = 2/t0, dda = 2/(t0*t0);
        }
        void evolution  ( double** f, double** df, int output_step ) override;
};


class LeapFrogSelf final: public LeapFrogBase
{
    private:
        void evolScale       ( double** f, const double h );
        void evolScaleDerivs ( const double h );

    public:
        LeapFrogSelf     ( std::shared_ptr<Model> model, int precision, double** f, double** df )
        : LeapFrogBase(model, precision)
        {
            double C = 0;
            for( int n = 0; n < num_fields; ++n ){
                #pragma omp parallel for reduction(+:C) schedule( static ) num_threads( num_threads )
                for( int i = 0; i < N; ++i ){
                    #if   DIMENSION == 1
                        int idx = i;
                        C += gradientEnergy(f[n], i) + 2*pow(a,4)*_model->V(f, n, idx, a);
                    #elif DIMENSION == 2
                        for( int j = 0; j < N; ++j ){
                            int idx = i*N+j;
                            C += gradientEnergy(f[n], i, j) + 2*pow(a,4)*_model->V(f, n, idx, a);
                        }
                    #elif DIMENSION == 3
                        for( int j = 0; j < N; ++j ){
                            for( int k = 0; k < N; ++k ){
                                int idx = (i*N+j)*N+k;
                                C += gradientEnergy(f[n], i, j, k) + 2*pow(a,4)*_model->V(f, n, idx, a);
                            }
                        }
                    #endif
                }
            }
            for( int d = 0; d < DIMENSION; ++d ) C /= N;
            da = sqrt(C/6);
        }

        void evolution       ( double** f, double** df, int output_step ) override;
};

#endif