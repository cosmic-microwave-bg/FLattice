#include <cmath>
#include "energy.hpp"



void Energy::calculate ( double**f, double** df )
{
    for( int n = 0; n < num_fields; ++n ){
        int i = 0, j = 0, k = 0;
        #pragma omp parallel for schedule( static ) num_threads ( num_threads )
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
                        if( expansion == Expansion::no_expansion )
                            _data[n][idx] = pow(df[n][idx], 2)/2
                                          + gradientEnergy(f[n], i, j, k)/2 + _model->V(f, n, idx);
                        else
                            _data[n][idx] = ( pow(df[n][idx] - f[n][idx]*da/a, 2)/2
                                          + gradientEnergy(f[n], i, j, k) )/(2*pow(a,4)) + _model->V(f, n, idx, a);
            #if DIMENSION == 3
                    }
            #endif
            #if DIMENSION >= 2
                }
            #endif
        }
    }
    
    #pragma omp parallel for schedule( static ) num_threads ( num_threads )
    for( int i = 0; i < N; ++i ){
        int idx = i;
        #if DIMENSION >= 2
        for( int j = 0; j < N; ++j ){
            idx = i*N+j;
        #endif
        #if DIMENSION == 3
            for( int k = 0; k < N; ++k ){
                idx = (i*N+j)*N+k;
        #endif
                for( int n = 0; n < num_fields; ++n ){
                    if( n == 0 ) _data_tot[idx]  = _data[n][idx];
                    else         _data_tot[idx] += _data[n][idx];
                }
        #if DIMENSION == 3
            }
        #endif
        #if DIMENSION >= 2
        }
        #endif
    }
}