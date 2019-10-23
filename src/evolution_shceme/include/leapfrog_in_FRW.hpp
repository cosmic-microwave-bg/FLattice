#ifndef _LEAPFROG_WITH_EXPANSION_H_
#define _LEAPFROG_WITH_EXPANSION_H_

#include <cmath>
#include "leapfrog_base.hpp"
#include "calculation.hpp"



class LeapFrogRadMat final: public LeapFrogBase
{
    private:
        double _index; // index of scale factor
        double _coeff; // coefficient of scale factor

    public:
        LeapFrogRadMat  ( int precision ): LeapFrogBase(precision) 
        {
            /** da is unnecessary when you use pi
            if(expansion == Expansion::rad_domination){
                if(A == 0) da = 1/(2*t0); // = c
                if(A == 1) da = c;
            }
            if(expansion == Expansion::mat_domination){
                if(A == 0) da = 2/(3*t0); // = c
                if(A == 1) da = pow(c, 4/3)*t0;
            }
            */
            if(A == 0){ // physical time
                _index = (expansion==Expansion::rad_domination)? 1./2: 2./3;
                _coeff = (expansion==Expansion::rad_domination)? 2*c: 3/2*c;
            }
            if(A == 1){ // conformal time
                _index = (expansion==Expansion::rad_domination)? 1: 2;
                _coeff = (expansion==Expansion::rad_domination)? c: pow(c, 2/3)/2;
            }
            t = t0 = 1/_coeff;
        }

        void evolution  ( double** f, double** df, int output_step ) override;
};


class LeapFrogSelf final: public LeapFrogBase
{
    private:
        void evolScale  ( double** f, const double h );
        void calcDa     ( double** f, double** df );

    public:
        LeapFrogSelf    ( int precision, double** f, double** df )
        : LeapFrogBase(precision) 
        {
            double rho = 0;
            // calculate energy density to determine da
            for( int n = 0; n < num_fields; ++n ){
                int i = 0, j = 0, k = 0;
                #pragma omp parallel for reduction(+:rho) schedule(static) num_threads (num_threads)
                for( i = 0; i < N; ++i ){
                    int idx = i;
                #if DIMENSION >= 2
                    for( j = 0; j < N; ++j ){
                        idx = i*N+j;
                #endif
                #if DIMENSION == 3
                        for( k = 0; k < N; ++k ){
                            idx = (i*N+j)*N+k;
                #endif
                            rho += pow(df[n][idx]/pow(a,D), 2)/2 + gradientEnergy(f[n], i, j, k)/pow(a,2) + V(f, n, idx);
                #if DIMENSION == 3
                        }
                #endif
                #if DIMENSION >= 2
                    }
                #endif
                }
            }
            rho /= pow(N, DIMENSION);

            // Set t0 from Hubble parameter
            da = pow(a, A+1)*R*sqrt(2*rho/(D*(D-1)));
            t = t0 = 1/da;
        } 

        void evolution  ( double** f, double** df, int output_step ) override;
};

#endif