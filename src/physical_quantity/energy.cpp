#include <cmath>
#include "energy.hpp"
#include "model.hpp"



void Energy::calculate (double**f, double** df)
{
    for( int n = 0; n < num_fields; ++n ){
        int i = 0, j = 0, k = 0;
        #pragma omp parallel for schedule(static) num_threads (num_threads)
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
                        _data[n][idx] = pow(df[n][idx]/pow(a,D), 2)/2 + gradientEnergy(f[n], i, j, k)/pow(a,2) + V(f, n, idx);
            #if DIMENSION == 3
                    }
            #endif
            #if DIMENSION >= 2
                }
            #endif
        }
    }
    
    #pragma omp parallel for schedule(static) num_threads (num_threads)
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
