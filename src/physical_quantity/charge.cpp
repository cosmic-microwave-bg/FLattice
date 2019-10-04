#include "charge.hpp"


void Charge::calculate(double** f, double** df)
{
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
                _data_tot[idx] = (f[0][idx]*df[1][idx] - df[0][idx]*f[1][idx]);
                if( expansion != Expansion::no_expansion ) _data_tot[idx] /= pow(a, D);
        #if DIMENSION == 3
            }
        #endif
        #if DIMENSION >= 2
        }
        #endif
    }
}
