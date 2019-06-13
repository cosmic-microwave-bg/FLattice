/**
 * @file calculation.hpp
 * @brief The functions calculating average, variance, gradient energy are defined.
 */
#ifndef _CALCULATION_H_
#define _CALCULATION_H_

#include <cmath>
#include "parameter.hpp"



template <typename dataType>
dataType calculateAverage( dataType* data, int N, int dimension )
{
    dataType average = 0;

    #pragma omp parallel for reduction(+:average) schedule( static ) num_threads ( num_threads )
    for( int i = 0; i < N; ++i ){
        switch( dimension )
        {
            case 1:
                int idx = i;
                average += data[idx];
                break;
            case 2:
                for( int j = 0; j < N; ++j ){
                    int idx = i*N+j;
                    average += data[idx];
                }
                break;
            case 3:
                for( int j = 0; j < N; ++j )
                    for( int k = 0; k < N; ++k ){
                        int idx = (i*N+j)*N+k;
                        average += data[idx];
                    }
                break;
        }
    }
    for( int i = 0; i < dimension; ++i ) average /= N;
    
    return average;
}


template <typename dataType>
dataType calculateVariance( dataType* data, dataType average, int N, int dimension )
{
    dataType variance = 0;

    #pragma omp parallel for reduction(+:variance) schedule( static ) num_threads ( num_threads )
    for( int i = 0; i < N; ++i ){
        switch( dimension )
        {
            case 1:
                int idx = i;
                variance += pow( data[idx]-average, 2 );
                break;
            case 2:
                for( int j = 0; j < N; ++j ){
                    int idx = i*N+j;
                    variance += pow( data[idx]-average, 2 );
                }
                break;
            case 3:
                for( int j = 0; j < N; ++j ){
                    for( int k = 0; k < N; ++k ){
                        int idx = (i*N+j)*N+k;
                        variance += pow( data[idx]-average, 2 );
                    }
                }
                break;
        }
    }
    for( int i = 0; i < dimension; ++i ) variance /= N;
    
    return variance;
}

template <typename dataType=double>
dataType gradientEnergy( double* f, int i, int j=0, int k=0 )
{
    dataType gradient_energy = 0;
    
    int ip1 = (i ==N-1) ?     0: i+1;
    int ip2 = (i >=N-2) ? i-N+2: i+2;
    int im1 = (i ==  0) ?   N-1: i-1;
    int im2 = (i <   2) ? i+N-2: i-2;

    int jp1 = (j == N-1)?     0: j+1;
    int jp2 = (j >= N-2)? j-N+2: j+2;
    int jm1 = (j ==   0)?   N-1: j-1;
    int jm2 = (j <    2)? j+N-2: j-2;

    #if DIMENSION == 1
        gradient_energy +=  pow( ( - f[ip2]  + 8*f[ip1]  - 8*f[im1] + f[im2] ) / (12*dx), 2 );
    #elif DIMENSION == 2
        gradient_energy += pow( ( - f[ip2*N+j] + 8*f[ip1*N+j] - 8*f[im1*N+j] + f[im2*N+j] ) / (12*dx), 2 );
        gradient_energy += pow( ( - f[i*N+jp2] + 8*f[i*N+jp1] - 8*f[i*N+jm1] + f[i*N+jm2] ) / (12*dx), 2 );
    #elif DIMENSION == 3
        int kp1 = (k == N-1)?     0: k+1;
        int kp2 = (k >= N-2)? k-N+2: k+2;
        int km1 = (k ==   0)?   N-1: k-1;
        int km2 = (k <    2)? k+N-2: k-2;                        
        gradient_energy += pow( ( - f[(ip2*N+j)*N+k] + 8*f[(ip1*N+j)*N+k] - 8*f[(im1*N+j)*N+k] + f[(im2*N+j)*N+k] ) / (12*dx), 2 );
        gradient_energy += pow( ( - f[(i*N+jp2)*N+k] + 8*f[(i*N+jp1)*N+k] - 8*f[(i*N+jm1)*N+k] + f[(i*N+jm2)*N+k] ) / (12*dx), 2 );
        gradient_energy += pow( ( - f[(i*N+j)*N+kp2] + 8*f[(i*N+j)*N+kp1] - 8*f[(i*N+j)*N+km1] + f[(i*N+j)*N+km2] ) / (12*dx), 2 );
    #endif

    return gradient_energy;
}

#endif