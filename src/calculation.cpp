/**
 * @file  calculation.cpp
 * @brief See the header file for details.
 */
#include <cmath>

#include "calculation.hpp"
#include "parameter.hpp"



void CalculateBase::calculate()
{
    for( int n = 0; n < num_fields; ++n ){
        #pragma omp parallel for schedule( static ) num_threads ( num_threads )
        for( int i = 0; i < N; ++i ){
            #if   DIMENSION == 1 
                calculateData(n, i);
            #elif DIMENSION == 2
                for( int j = 0; j < N; ++j ) calculateData(n, i, j);
            #elif DIMENSION == 3
                for( int j = 0; j < N; ++j ){
                    for( int k = 0; k < N; ++k ) calculateData(n, i, j, k);
                }
            #endif
        }
    }
}


void CalculateEnergy::calculateData( int n, int i, int j, int k )
{
    double a  = LeapFrogExpansion::_a;
    double da = LeapFrogExpansion::_da;

    int idx = (k*N+j)*N+i;
    if( expansion == Expansion::no_expansion ) _data[n][idx] = pow(_df[n][idx], 2)/2 + gradientEnergy(_f[n], i, j, k)          + _model->V(_f, n, idx);
    else          _data[n][idx] = pow(_df[n][idx]*a - _f[n][idx]*da, 2)/(2*pow(a,4)) + gradientEnergy(_f[n], i, j, k)/pow(a,4) + _model->aV(_f, a, n, idx);
    
}


double CalculateEnergy::gradientEnergy( const double* f, int i, int j, int k ) const
{
    double gradient_energy = 0;
    
    int ip1 = (i == N-1)?     0: i+1;
    int ip2 = (i >= N-2)? i-N+2: i+2;
    #ifdef SPHERICAL_SYM
    int im1 = abs( i-1 );
    int im2 = abs( i-2 );
    #else
    int im1 = (i == 0)?   N-1: i-1;
    int im2 = (i <  2)? i+N-2: i-2;
    #endif

    #if DIMENSION == 1
        #ifdef SPHERICAL_SYM
        int idx = i;
        if( idx == 0 ){
            grad[0] = 0;
        }else if( idx == N-2 ){
            grad[0] = ( f[ip1] - f[im1] ) / (2*dx);
        }else if( idx == N-1 ){
            grad[0] = ( f[im2] - 4*f[im2] + 3*f[idx]  ) / (2*dx);
        }else{
            grad[0] = ( - f[ip2]  + 8*f[ip1]  - 8*f[im1] + f[im2] ) / (12*dx);
        }
        #else
        gradient_energy +=  pow( ( - f[ip2]  + 8*f[ip1]  - 8*f[im1] + f[im2] ) / (12*dx), 2 );
        #endif
    #elif DIMENSION == 2
            int jp1 = (j == N-1)?     0: j+1;
            int jp2 = (j >= N-2)? j-N+2: j+2;
            int jm1 = (j ==   0)?   N-1: j-1;
            int jm2 = (j <    2)? j+N-2: j-2;
            gradient_energy += pow( (- f[ip2*N+j] + 8*f[ip1*N+j] - 8*f[im1*N+j] + f[im2*N+j]) / (12*dx), 2 );
            gradient_energy += pow( (- f[i*N+jp2] + 8*f[i*N+jp1] - 8*f[i*N+jm1] + f[i*N+jm2]) / (12*dx), 2 );
    #elif DIMENSION == 3
            int jp1 = (j == N-1)?     0: j+1;
            int jp2 = (j >= N-2)? j-N+2: j+2;
            int jm1 = (j ==   0)?   N-1: j-1;
            int jm2 = (j <    2)? j+N-2: j-2;
                int kp1 = (k == N-1)?     0: k+1;
                int kp2 = (k >= N-2)? k-N+2: k+2;
                int km1 = (k ==   0)?   N-1: k-1;
                int km2 = (k <    2)? k+N-2: k-2;                        
                gradient_energy += pow( (- f[(ip2*N+j)*N+k] + 8*f[(ip1*N+j)*N+k] - 8*f[(im1*N+j)*N+k] + f[(im2*N+j)*N+k]) / (12*dx), 2 );
                gradient_energy += pow( (- f[(i*N+jp2)*N+k] + 8*f[(i*N+jp1)*N+k] - 8*f[(i*N+jm1)*N+k] + f[(i*N+jm2)*N+k]) / (12*dx), 2 );
                gradient_energy += pow( (- f[(i*N+j)*N+kp2] + 8*f[(i*N+j)*N+kp1] - 8*f[(i*N+j)*N+km1] + f[(i*N+j)*N+km2]) / (12*dx), 2 );
    #endif

	return gradient_energy / 2;
}